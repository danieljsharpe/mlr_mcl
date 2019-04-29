/*
C++ code for multi-level regularised Markov clustering (MLR-MCL) of a kinetic transition network
Need two files: "ts_conns.dat" (two-column format: min1, min2) and "ts_weights.dat" (one-column format)
in the current directory

Compile with:
g++ -std=c++11 mlr_mcl.cpp ktn.cpp read_ktn.cpp utils.cpp -I /usr/include/python2.7/ -L /usr/include/python2.7/Python.h -lpython2.7 -o mlr_mcl

Execute with e.g.:
./mlr_mcl 22659 34145 1.5 0.5 10 4 1.E-06 1.E-06 19 2000

Daniel J. Sharpe
April 2019
*/

#include "read_ktn.h"
#include "ktn.h"
#include "utils.h"
#include <iostream>
#include <array>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>
#include <map>
#include <stdlib.h>
#include <boost/python.hpp>
#include <Python.h>

using namespace std;

// class for multi-level regularised Markov clustering (MLR-MCL)
class MLR_MCL {

    typedef vector<vector<pair<int,int>>> nodemap_vec;
    typedef pair<vector<pair<double,int>>,vector<int>> Csr_mtx; // matrix in CSR sparse format

    public:
    MLR_MCL(double,double,int,int,double,double,int,int);
    ~MLR_MCL();
    void run_mcl(Network&);

    double r; double b; double eps; double tau;
    int n_C; int n_cur; int seed; int min_C;

    private:
    void coarsen_graph(Network&);
    void heavy_edge_matching(Network&,int);
    void get_matrix_exp(Network&);
    Csr_mtx get_init_sparse_mtx(Network&);
    vector<int> get_col_idcs(Network&,vector<int>,int);
    void prune();
    void regularise();
    void expand();
    void inflate();
    void project_flow();
    void interpret_clust();

    vector<int> g_C_size; // size of coarsened graph at each step
    nodemap_vec nodemap; // mapping of nodes to higher levels in multi-level coarsening procedure
};

MLR_MCL::MLR_MCL(double d1,double d2,int i1,int i2,double d3, double d4, int i3, int i4) : \
                r(d1),b(d2),n_C(i1),n_cur(i2),eps(d3),tau(d4),seed(i3),min_C(i4) {

    nodemap.resize(n_C);
};

MLR_MCL::~MLR_MCL() {};

void MLR_MCL::run_mcl(Network &ktn) {
    coarsen_graph(ktn);
    get_matrix_exp(ktn);
}

// heavy edge matching for graph coarsening
void MLR_MCL::heavy_edge_matching(Network &ktn, int i_C) {
    Edge *edgeptr;
    vector<int> node_ids(ktn.tot_nodes);
    iota(begin(node_ids),end(node_ids),1);
    random_shuffle(begin(node_ids),end(node_ids));
//    vector<int> node_ids;
//    if (ktn.n_nodes==8) { node_ids = {1,3,5,7,2,4,6,8}; } // quack TEST ONLY
//    else if (ktn.n_nodes==4) { node_ids = {1,5,3,7,1,1,1,} ; }
    int matchnode_id;
    for (int i=0;i<ktn.tot_nodes;i++) {
        ktn.min_nodes[i].hem_flag = false; }
    for (int i=0;i<ktn.tot_nodes;i++) {
//        cout << "\n iteration no. " << i+1 << " chosen node... " << ktn.min_nodes[node_ids[i]-1].min_id << endl;
        // check if node has already been matched or deleted
        if ((ktn.min_nodes[node_ids[i]-1].hem_flag) || (ktn.min_nodes[node_ids[i]-1].deleted)) { continue; }
        else { ktn.min_nodes[node_ids[i]-1].hem_flag = true; }
//        cout << "  has not been deleted or matched..." << endl;
        // find the neighbour FROM node i with largest weight
        edgeptr = ktn.min_nodes[node_ids[i]-1].top_from;
        double best_w = -numeric_limits<double>::infinity(); matchnode_id = 0;
        if (edgeptr != nullptr) {
        do {
            if ((edgeptr->to_node->hem_flag) || \
                (ktn.min_nodes[edgeptr->to_node->min_id-1].deleted)) {
                edgeptr = edgeptr->next_from; continue; }
            if (edgeptr->w > best_w) {
                best_w = edgeptr->w;
                matchnode_id = edgeptr->to_node->min_id;
            }
            edgeptr = edgeptr->next_from;
        } while (edgeptr!=nullptr);
        }
        if (matchnode_id != 0) { // found a matching
//            cout << "    highest weight FROM nbr is node " << matchnode_id << " w " << best_w << endl;
            ktn.min_nodes[matchnode_id-1].hem_flag = true;
            ktn.merge_nodes(ktn.min_nodes[node_ids[i]-1].min_id-1,ktn.min_nodes[matchnode_id-1].min_id-1);
//            cout << "    finished merging nodes" << endl;
            // update nodemaps
            nodemap[i_C].emplace_back(make_pair(node_ids[i],matchnode_id));
        } else { ktn.del_node(node_ids[i]-1); }
    }
}

// function for coarsening of the graph
void MLR_MCL::coarsen_graph(Network &ktn) {
    int i=0;
    if (ktn.n_nodes > min_C) {
    do {
        heavy_edge_matching(ktn,i);
        i++;
        cout << ">>>> n_nodes IS NOW: " << ktn.n_nodes << endl;
        g_C_size.emplace_back(ktn.n_nodes);
    } while ((i < n_C) && (ktn.n_nodes > min_C));
    }
    cout << "finished graph coarsening after " << i << " iterations" << endl;
}

// call external Python script to compute the matrix exponential of the coarsened graph.
// Return the initial, coarsened column-stochastic flow matrix
void MLR_MCL::get_matrix_exp(Network &ktn) {
/*
    Edge *edgeptr;
    cout << "looping over all FROM neighbours: " << endl;
    for (int i=0;i<ktn.tot_nodes;i++) {
        if (ktn.min_nodes[i].deleted) { continue; }
        cout << " node " << ktn.min_nodes[i].min_id << endl;
        edgeptr = ktn.min_nodes[i].top_from;
        if (edgeptr!=nullptr) {
        do {
            cout << "  FROM " << edgeptr->from_node->min_id << " TO " << edgeptr->to_node->min_id << " w " << \
                    edgeptr->w << "  TO node is deleted? " << edgeptr->to_node->deleted << endl;
            edgeptr = edgeptr->next_from; } while (edgeptr!=nullptr);
        } else { cout << "  node has no FROM nbrs!" << endl; }
    }
*/
    Csr_mtx k_mtx_sp = get_init_sparse_mtx(ktn); // sparse representation of transition rate matrix
    cout << "Calculated sparse transition rate matrix" << endl;
    int ret_val;
    const char *py_script_name = "calc_matrix_exp"; const char *py_func_name = "calc_expm";
    setenv("PYTHONPATH",".",1); // set environment variable correctly (here, calc_matrix_exp.py is in current dir)
    PyObject *pName=NULL, *pModule=NULL, *pFunc=NULL; // names of Python script, module and function, respectively
    PyObject *pArgs=NULL, *pRetval=NULL; // names of arguments for and return value of the Python function, respectively
    // PyObject objs for defining arguments to the Python function
    PyObject *pValue=NULL, *pValue2=NULL, *pValue3=NULL, *pValue4=NULL;
    Py_Initialize(); // initialise the interpreter
    pName = PyString_FromString(py_script_name); // name of Python script
    if (pName==NULL) { goto error; }
    pModule = PyImport_Import(pName); // import corresponding module
    Py_DECREF(pName);
    if (pModule==NULL) { goto error; }
    pFunc = PyObject_GetAttrString(pModule,py_func_name); // retrieve desired object from module
    if (pFunc==NULL) { goto error; }
    if (pFunc && PyCallable_Check(pFunc)) { // Python object exists and is callable
        pArgs = PyTuple_New(4); // tuple of args. Here, pFunc takes four args; 3 arrays and 1 double
        pRetval = PyInt_FromLong(0); // dummy
        // each arg is a Python object. Here, we want to pass 3 lists describing the sparse matrix to the Python function object
        pValue = PyList_New(k_mtx_sp.first.size());
        pValue2 = PyList_New(k_mtx_sp.first.size());
        pValue3 = PyList_New(k_mtx_sp.second.size());
        for (int i=0;i<k_mtx_sp.first.size();i++) {
            PyList_SetItem(pValue,i,PyFloat_FromDouble(k_mtx_sp.first[i].first));
            PyList_SetItem(pValue2,i,PyInt_FromLong(k_mtx_sp.first[i].second));
        }
        for (int i=0;i<k_mtx_sp.second.size();i++) {
            PyList_SetItem(pValue3,i,PyInt_FromLong(k_mtx_sp.second[i]));
        }
        pValue4 = PyFloat_FromDouble(tau);
        // set values of the pArgs tuple
        PyTuple_SetItem(pArgs,0,pValue); PyTuple_SetItem(pArgs,1,pValue2); PyTuple_SetItem(pArgs,2,pValue3);
        PyTuple_SetItem(pArgs,3,pValue4);
        pRetval = PyObject_CallObject(pFunc,pArgs); // call the Python function
        cout << "I'm back outside of the Python script" << endl;
        if (pValue==NULL) { goto error; } // Python function has a return value, so should not return NULL
        ret_val = int(PyFloat_AsDouble(pRetval)); // extract the return value of the Python function
        Py_DECREF(pValue); Py_DECREF(pValue2); Py_DECREF(pValue3); Py_DECREF(pValue4);
        Py_DECREF(pFunc); Py_DECREF(pModule); Py_DECREF(pArgs); Py_DECREF(pRetval);
    } else { cout << "Error: could not find Python function object" << endl; exit(EXIT_FAILURE); }
    error: // check if error has occurred in use of Python interpreter
        if (PyErr_Occurred()) {
            cout << "Fatal error in Python interpreter" << endl; exit(EXIT_FAILURE); }
    Py_Finalize(); // finished with Python interpreter
    cout << "Return value from py function was: " << ret_val << endl;
}

// construct an initial sparse matrix representation of the coarsened network
MLR_MCL::Csr_mtx MLR_MCL::get_init_sparse_mtx(Network &ktn) {

    // elements of the transition rate matrix, row-major order, and corresponding column indices
    vector<pair<double,int>> k_elems_cols;
    vector<int> k_rl; // cumulative row lengths of transition rate matrix
    vector<double> k_elems; vector<int> k_minids; // temporary vectors for matrix elems and min_id's
    map<int,int> k_min2col; // mapping of min_id's to column indices
    Edge *edgeptr;
    // set the elements of the (flattened) transition rate matrix
    int z=0, dconn=0;
    for (int i=0;i<ktn.tot_nodes;i++) {
        int rl=0;
        if (ktn.min_nodes[i].deleted) { continue; }
        edgeptr = ktn.min_nodes[i].top_to;
        bool edge_exist = false;
        if (edgeptr!=nullptr) {
        do {
            // note that the column indices/min id's are not initially in order
            if (edgeptr->from_node->deleted) { edgeptr = edgeptr->next_to; continue; }
            if (!edge_exist) { edge_exist = true; }
            k_elems.emplace_back(edgeptr->w); k_minids.emplace_back(edgeptr->from_node->min_id);
            try {
                int dummy = k_min2col.at(edgeptr->from_node->min_id);
            } catch (const out_of_range& oor) {
                k_min2col[edgeptr->from_node->min_id] = -1; z++; } // initial dummy value in map
            rl++;
            edgeptr = edgeptr->next_to;
        } while (edgeptr != nullptr);
        }
        if (!edge_exist) { dconn++; } // this node is not deleted but is disconnected
        k_rl.emplace_back(rl);
    }
    cout << "number of unique col_idcs: " << z << endl;
    cout << "number of disconnected nodes: " << dconn << endl;
//    for (int i=0;i<k_minids.size();i++) { cout << k_minids[i] << endl; }
    for (vector<int>::iterator it=k_rl.begin()+1;it!=k_rl.end();it++) {
        *it += *(it-1); }
    vector<int> k_minids_sort = k_minids;
    sort(k_minids_sort.begin(),k_minids_sort.end());
    int i=0;
    for (auto minid: k_minids_sort) {
        if (k_min2col[minid]==-1) { k_min2col[minid] = i; i++; } }
    for (int i=0;i<k_minids.size();i++) {
        k_elems_cols.emplace_back(make_pair(k_elems[i],k_min2col[k_minids[i]])); }
    Csr_mtx sparse_mtx = make_pair(k_elems_cols,k_rl);
    return sparse_mtx;
}

/*
// given the min IDs of a nodemap, get the corresponding column indices of the transition matrix
vector<pair<,int,int>> MLR_MCL::get_col_idcs(Network& ktn, int n_it) {
    vector<pair<int,int>> nodemap_col(ktn.nodemap_vec[n_it].size());
    vector<pair<int,int>> nodemap_copy = ktn.nodemap_vec[n_it];
    sort(nodemap_copy.begin(),nodemap_copy.end(),[](pair<int,int> pair1, pair<int,int> pair2) \
         { return (pair1.first < pair2.first); }
    
    
    return nodemap_col;
}
*/

// prune small values from the transition matrix
void MLR_MCL::prune() {

}

// regularisation operation for the transition matrix
void MLR_MCL::regularise() {

}

// expansion operation for the transition matrix
void MLR_MCL::expand() {

}

// inflation operation for the transition matrix
void MLR_MCL::inflate() {

}

// use nodemaps to undo the coarsening of the graph
void MLR_MCL::project_flow() {

}

// after the clustering procedure has finished, interpret the flow matrix as a clustering
// characterised by attractors
void MLR_MCL::interpret_clust() {

}

int main(int argc, char** argv) {

    // initialise
    int nmin = stoi(argv[1]); // no. of minima
    int nts = stoi(argv[2]); // no. of transition states
    double r = stod(argv[3]); // inflation parameter
    double b = stod(argv[4]); // balance parameter
    int n_C = stoi(argv[5]); // max no. of coarsenings
    int n_cur = stoi(argv[6]); // no. of curtailed MCL iterations for coarsened graphs
    double eps = stod(argv[7]); // threshold for pruning
    double tau = stod(argv[8]); // lag time for estimating transition matrix from transition rate matrix
    int seed = stoi(argv[9]); // random seed
    int min_C; // min. no. of nodes in coarsened graph
    if (argc > 10) { min_C = stoi(argv[10]); } else { min_C = 0; }

    vector<pair<int,int>> ts_conns = Read_ktn::read_ts_conns(nts);
    vector<double> ts_weights = Read_ktn::read_ts_weights(nts);

    Network ktn(nmin,nts);
    Network::setup_network(ktn,nmin,nts,ts_conns,ts_weights);
    ts_conns.resize(0); ts_weights.resize(0);

//    run_debug_tests(ktn);

    MLR_MCL mcl_obj (r,b,n_C,n_cur,eps,tau,seed,min_C);
    mcl_obj.run_mcl(ktn);

    return 0;
}
