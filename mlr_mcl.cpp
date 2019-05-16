/*
C++ code for multi-level regularised Markov clustering (MLR-MCL) of a kinetic transition network
Need two files: "ts_conns.dat" (two-column format: min1, min2) and "ts_weights.dat" (one-column format)
in the current directory

Compile with:
g++ -std=c++11 mlr_mcl.cpp ktn.cpp read_ktn.cpp utils.cpp -I /usr/include/python2.7/ -L /usr/include/python2.7/Python.h -lpython2.7 -o mlr_mcl

Execute with e.g.:
./mlr_mcl 22659 34145 1.5 0.5 10 4 100 1.E-12 1.E-10 19 1000

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
#include <iterator>
#include <stdlib.h>
#include <boost/python.hpp>
#include <Python.h>

using namespace std;

/* class for multi-level regularised Markov clustering (MLR-MCL) */
class MLR_MCL {

    typedef vector<map<int,int>> nodemap_vec;
    typedef pair<vector<pair<double,int>>,vector<int>> Csr_mtx; // matrix in CSR (or CSC) sparse format

    public:
    MLR_MCL(double,double,int,int,int,double,double,int,int);
    ~MLR_MCL();
    void run_mcl(Network&);

    double r; double b; double eps; double tau;
    int n_C; int n_cur; int max_it; int seed; int min_C;

    private:
    void mcl_main_ops(Csr_mtx&,const Csr_mtx&);
    void coarsen_graph(Network&);
    void heavy_edge_matching(Network&,int);
    Csr_mtx get_matrix_exp(Network&);
    Csr_mtx get_init_sparse_mtx(Network&);
    void regularise(Csr_mtx&,const Csr_mtx&);
    void inflate(Csr_mtx&);
    void prune_renormalise(Csr_mtx&);
    Csr_mtx get_reg_mtx(const Csr_mtx&);
    Csr_mtx project_flow(Network&,Csr_mtx&,map<int,int>);
    void interpret_clust(const Csr_mtx&);

    vector<int> g_C_size; // size of coarsened graph at each step
    nodemap_vec nodemap; // mapping of nodes to higher levels in multi-level coarsening procedure (pairwise merging)
    map<int,int> idxmap; // mapping of min ID's of coarsened graph to indices of corresponding transition matrix
    vector<int> idxlist; // ordered list of min ID's (posns in list correspond to indices of transition matrix)
};

MLR_MCL::MLR_MCL(double d1,double d2,int i1,int i2,int i3,double d3,double d4,int i4,int i5) : \
                r(d1),b(d2),n_C(i1),n_cur(i2),max_it(i3),eps(d3),tau(d4),seed(i4),min_C(i5) {

    nodemap.resize(n_C);
};

MLR_MCL::~MLR_MCL() {};

/* main loop to drive multi-level regularised Markov clustering */
void MLR_MCL::run_mcl(Network &ktn) {
    coarsen_graph(ktn);
    Csr_mtx t_mtx_sp = get_matrix_exp(ktn); // transition matrix (CSR format)
    Csr_mtx tG_mtx_sp = get_reg_mtx(t_mtx_sp); // regularisation matrix (CSC format)

//    cout << "test printing first elems of initial transition mtx..." << endl;
//    for (int i=0;i<10;i++) {
//        cout << t_mtx_sp.first[i].first << "  " << t_mtx_sp.first[i].second << "  " << t_mtx_sp.second[i] << endl; }

    // run curtailed MLR-MCL
    for (int i=g_C_size.size()-1;i>=1;i--) {
        cout << ">>>>> running curtailed MLR-MCL on coarsened graph at level " << i+1 << endl;
/*
        cout << "  length of nodemap: " << nodemap[i].size() << endl;
        map<int,int>::iterator it_map;
        int k=0;
        for (it_map=nodemap[i].begin();it_map!=nodemap[i].end();it_map++) {
            cout << "  node 1: " << it_map->first << " maps to node 2: " << it_map->second << endl;
            k++; if (k>10) { break; }
        }
*/
        for (int j=0;j<n_cur;j++) { mcl_main_ops(t_mtx_sp,tG_mtx_sp); }
        t_mtx_sp = project_flow(ktn,t_mtx_sp,nodemap[i]); // refined transition matrix
        cout << "returned projected t_mtx_sp" << endl;
        tG_mtx_sp = get_reg_mtx(t_mtx_sp);
    }
//    for (int i=0;i<max_it;i++) { mcl_main_ops(t_mtx_sp,tG_mtx_sp); }
    interpret_clust(t_mtx_sp);
}

/* operations of Markov clustering main loop */
void MLR_MCL::mcl_main_ops(Csr_mtx &t_mtx_sp, const Csr_mtx &tG_mtx_sp) {
    regularise(t_mtx_sp,tG_mtx_sp);
    inflate(t_mtx_sp);
    cout << "no. elems before prune " << t_mtx_sp.first.size() << endl;
    prune_renormalise(t_mtx_sp);
    cout << "no. elems after prune " << t_mtx_sp.first.size() << endl;
}

/* heavy edge matching for graph coarsening */
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
            // update nodemap
            nodemap[i_C][node_ids[i]] = matchnode_id;
        } else { ktn.del_node(node_ids[i]-1); }
    }
}

/* function for coarsening of the graph */
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

/* call external Python script to compute the matrix exponential of the coarsened graph.
   Return the initial, coarsened column-stochastic flow matrix in CSR sparse format */
MLR_MCL::Csr_mtx MLR_MCL::get_matrix_exp(Network &ktn) {
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
    cout << "\n\n\n\n" << endl;
*/
    cout << "Calculating sparse transition rate matrix" << endl;
    Csr_mtx k_mtx_sp = get_init_sparse_mtx(ktn); // sparse representation of transition rate matrix
    cout << "Setting up Python interpreter" << endl;
    vector<pair<double,int>> Tspci; vector<int> Trl;
    const char *py_script_name = "calc_matrix_exp"; const char *py_func_name = "calc_expm";
    bool err_flag = false;
    setenv("PYTHONPATH",".",1); // set environment variable correctly (here, calc_matrix_exp.py is in current dir)
    PyObject *pName=NULL, *pModule=NULL, *pFunc=NULL; // names of Python script, module and function, respectively
    PyObject *pArgs=NULL, *pRetval=NULL; // names of arguments for and return value of the Python function, respectively
    // PyObject objs for defining arguments to the Python function
    PyObject *pValue=NULL, *pValue2=NULL, *pValue3=NULL, *pValue4=NULL, *pValue5=NULL;
    PyObject *pList1=NULL, *pList2=NULL, *pList3=NULL; // used to process returned lists
    Py_Initialize(); // initialise the interpreter
    pName = PyString_FromString(py_script_name); // name of Python script
    if (pName==NULL) { err_flag = true; goto error; }
    pModule = PyImport_Import(pName); // import corresponding module
    Py_DECREF(pName);
    if (pModule==NULL) { err_flag = true; goto error; }
    pFunc = PyObject_GetAttrString(pModule,py_func_name); // retrieve desired object from module
    if (pFunc==NULL) { err_flag = true; goto error; }
    if (pFunc && PyCallable_Check(pFunc)) { // Python object exists and is callable
        pArgs = PyTuple_New(5); // tuple of args. Here, pFunc takes four args; 3 arrays and 1 double
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
        pValue5 = PyFloat_FromDouble(eps);
        // set values of the pArgs tuple
        PyTuple_SetItem(pArgs,0,pValue); PyTuple_SetItem(pArgs,1,pValue2); PyTuple_SetItem(pArgs,2,pValue3);
        PyTuple_SetItem(pArgs,3,pValue4); PyTuple_SetItem(pArgs,4,pValue5);
        pRetval = PyObject_CallObject(pFunc,pArgs); // call the Python function
        cout << "Returned transition matrix from Python script" << endl;
        // Python function has a return value, so should not return NULL, and should return a tuple of lists
        if (pRetval==NULL || !PyTuple_Check(pRetval)) { err_flag = true; goto error; }
        pList1 = PyTuple_GetItem(pRetval,0); pList2 = PyTuple_GetItem(pRetval,1); pList3 = PyTuple_GetItem(pRetval,2);
        if (!PyList_Check(pList1) || !PyList_Check(pList2) || !PyList_Check(pList3)) { err_flag = true; goto error; }
        // extract the return values of the Python function
        cout << "Transition matrix has " << PyList_Size(pList1) << " nonzero elems and dimension " << \
                PyList_Size(pList3) << endl;
        Tspci.resize(PyList_Size(pList1)); Trl.resize(PyList_Size(pList3));
        for (int i=0;i<PyList_Size(pList1);i++) {
            double Tsp_i = PyFloat_AsDouble(PyList_GetItem(pList1,i));
            int Tci_i = int(PyInt_AsLong(PyList_GetItem(pList2,i)));
            Tspci[i] = make_pair(Tsp_i,Tci_i); }
        for (int i=0;i<PyList_Size(pList3);i++) {
            Trl[i] = int(PyInt_AsLong(PyList_GetItem(pList3,i))); }
        Py_DECREF(pValue); Py_DECREF(pValue2); Py_DECREF(pValue3); Py_DECREF(pValue4);
        Py_DECREF(pFunc); Py_DECREF(pModule); Py_DECREF(pArgs); Py_DECREF(pRetval);
        Py_DECREF(pList1); Py_DECREF(pList2); Py_DECREF(pList3);
    } else { cout << "Error: could not find Python function object" << endl; exit(EXIT_FAILURE); }
    error: // check if error has occurred in use of Python interpreter
        if (PyErr_Occurred()) {
            cout << "Fatal error in Python interpreter" << endl; exit(EXIT_FAILURE);
        } else if (err_flag) {
            cout << "Fatal: manual error checking has raised a flag" << endl; exit(EXIT_FAILURE); }
    Py_Finalize(); // finished with Python interpreter
    Csr_mtx t_mtx_sp = make_pair(Tspci,Trl);
    return t_mtx_sp;
}

/* construct an initial sparse matrix representation of the coarsened network */
MLR_MCL::Csr_mtx MLR_MCL::get_init_sparse_mtx(Network &ktn) {

    // elements of the transition rate matrix, row-major order, and corresponding column indices
    vector<pair<double,int>> k_elems_cols;
    vector<int> k_rl; // cumulative row lengths of transition rate matrix
    vector<double> k_elems; vector<int> k_minids; // temporary vectors for matrix elems and column indices
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
            if (edgeptr->from_node->min_id-1==i) { throw Network::Ktn_exception(); }
            // note that the column indices/min id's are not initially in order
            if (edgeptr->from_node->deleted) { edgeptr = edgeptr->next_to; continue; }
            if (!edge_exist) { edge_exist = true; }
            k_elems.emplace_back(edgeptr->w); k_minids.emplace_back(edgeptr->from_node->min_id);
            try {
                int dummy = idxmap.at(edgeptr->from_node->min_id);
            } catch (const out_of_range& oor) {
                idxmap[edgeptr->from_node->min_id] = -1; z++; } // initial dummy value in map
            rl++;
            edgeptr = edgeptr->next_to;
        } while (edgeptr != nullptr);
        }
        if (!edge_exist) { ktn.min_nodes[i].deleted = true; dconn++; } // this node is not deleted but is disconnected
        else { k_rl.emplace_back(rl); }
    }
    cout << "number of unique col_idcs: " << z << endl;
    cout << "number of disconnected nodes: " << dconn << endl;
    for (vector<int>::iterator it=k_rl.begin()+1;it!=k_rl.end();it++) {
        *it += *(it-1); }
    idxlist.reserve(idxmap.size());
    transform(begin(idxmap),end(idxmap),back_inserter(idxlist), \
              [](decltype(idxmap)::value_type const& pair) { return pair.first; });
    int i=0;
    for (auto minid: idxlist) {
        idxmap[minid] = i; i++; }
    for (int i=0;i<k_elems.size();i++) {
        k_elems_cols.emplace_back(make_pair(k_elems[i],idxmap[k_minids[i]])); }
    Csr_mtx sparse_mtx = make_pair(k_elems_cols,k_rl);
    return sparse_mtx;
}

/* construct regularisation matrix and perform regularisation operation for the transition matrix */
void MLR_MCL::regularise(Csr_mtx &T_csr, const Csr_mtx &TG_csr) {

}

/* inflation operation for the transition matrix */
void MLR_MCL::inflate(Csr_mtx &T_csr) {
    for (int i=0;i<T_csr.first.size();i++) {
        T_csr.first[i].first = pow(T_csr.first[i].first,r); }
}

/* prune small values and renormalise columns of sparse transition matrix */
void MLR_MCL::prune_renormalise(Csr_mtx &T_csr) {

    vector<double> cum_sum(T_csr.second.size(),0.);
    int k=0, rn=0;
    vector<pair<double,int>>::iterator it_vec = T_csr.first.begin();
    while (it_vec!=T_csr.first.end()) { // prune and accumulate column sums
        if (it_vec->first < eps) {
            it_vec=T_csr.first.erase(it_vec);
            T_csr.second[rn]--;
        }
        else {
            cum_sum[it_vec->second] += it_vec->first; it_vec++; }
        if (k==T_csr.second[rn]) { rn++; }; k++;
    }
    k=0; rn=0;
    for (it_vec=T_csr.first.begin();it_vec!=T_csr.first.end();it_vec++) { // renormalise columns
        it_vec->first *= (1./cum_sum[it_vec->second]); }
}

/* given a transition matrix in CSR format, return a regularisation matrix in CSC format */
MLR_MCL::Csr_mtx MLR_MCL::get_reg_mtx(const Csr_mtx& T_csr) {

    Csr_mtx TG_csr;
    return TG_csr;
}

/* use nodemaps to undo the coarsening of the graph (refinement) */
MLR_MCL::Csr_mtx MLR_MCL::project_flow(Network &ktn, Csr_mtx &T_csr, map<int,int> curr_nodemap) {

    vector<int> idxlist_old = idxlist;
    map<int,int>::iterator it_map; vector<int>::iterator it_find;
    for (it_map=curr_nodemap.begin();it_map!=curr_nodemap.end();it_map++) { // update the list of indices
        if (!ktn.min_nodes[it_map->first-1].deleted) { // add second node of pair to list of indices
            idxlist.emplace_back(it_map->second); }
    }
    sort(idxlist.begin(),idxlist.end());
    int k=0;
    for (auto minid: idxlist) { // update the map of indices
        idxmap[minid] = k; k++;
    }
    // elements and columns of refined transition matrix in adjacency list format
    vector<vector<pair<double,int>>> T_ec_new(idxlist.size());
    vector<int> T_rl(idxlist.size(),0);
    k=0; int rn=0; // row no.
    vector<pair<double,int>>::iterator it_vec;
    // only two of the four corresponding elements in the refined matrix are non-zero
    for (it_vec=T_csr.first.begin();it_vec!=T_csr.first.end();it_vec++) {
        T_ec_new[idxmap[idxlist_old[rn]]].emplace_back(make_pair(it_vec->first, \
                 idxmap[idxlist_old[it_vec->second]]));
        T_ec_new[idxmap[idxlist_old[rn]]].emplace_back(make_pair(it_vec->first, \
                 curr_nodemap[idxmap[idxlist_old[it_vec->second]]]));
        if (k==T_csr.second[rn]) { rn++; }; k++;
    }
    k=0;
    for (auto rowvec: T_ec_new) {
        T_rl[k] = rowvec.size(); k++; }
    // flatten the refined transition matrix
    vector<pair<double,int>> T_ec = flatten<pair<double,int>>(T_ec_new);
    Csr_mtx T_csr_new = make_pair(T_ec,T_rl);
    return T_csr_new;
}

/* after the clustering procedure has finished, interpret the flow matrix as a clustering
   characterised by attractors */
void MLR_MCL::interpret_clust(const Csr_mtx &T_csr) {

}

int main(int argc, char** argv) {

    // initialise
    int nmin = stoi(argv[1]); // no. of minima
    int nts = stoi(argv[2]); // no. of transition states
    double r = stod(argv[3]); // inflation parameter
    double b = stod(argv[4]); // balance parameter
    int n_C = stoi(argv[5]); // max no. of coarsenings
    int n_cur = stoi(argv[6]); // no. of curtailed MCL iterations for coarsened graphs
    int max_it = stoi(argv[7]); // max. no. of iterations of R-MCL on the full graph
    double eps = stod(argv[8]); // threshold for pruning
    double tau = stod(argv[9]); // lag time for estimating transition matrix from transition rate matrix
    int seed = stoi(argv[10]); // random seed
    int min_C; // min. no. of nodes in coarsened graph (optional)
    int debug_flag; // run debug tests Y/N (optional)
    if (argc > 11) { min_C = stoi(argv[11]); } else { min_C = 0; }
    if (argc > 12) { debug_flag = stoi(argv[12]); } else { debug_flag = 0; }

    vector<pair<int,int>> ts_conns = Read_ktn::read_ts_conns(nts);
    vector<double> ts_weights = Read_ktn::read_ts_weights(nts);

    Network ktn(nmin,nts);
    Network::setup_network(ktn,nmin,nts,ts_conns,ts_weights);
    ts_conns.resize(0); ts_weights.resize(0);

    if (debug_flag) { run_debug_tests(ktn); exit(0); }

    MLR_MCL mcl_obj (r,b,n_C,n_cur,max_it,eps,tau,seed,min_C);
    mcl_obj.run_mcl(ktn);

    return 0;
}
