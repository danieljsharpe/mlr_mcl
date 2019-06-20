/*
C++ code for multi-level regularised Markov clustering (MLR-MCL) of a kinetic transition network
Need two files: "ts_conns.dat" (two-column format: min1, min2) and "ts_weights.dat" (one-column format)
in the current directory

Compile with:
g++ -std=c++11 mlr_mcl.cpp ktn.cpp utils.cpp quality.cpp -I /usr/include/python2.7/ -L /usr/include/python2.7/Python.h -lpython2.7 -o mlr_mcl

Execute with e.g.:
./mlr_mcl 22659 34145 1.2 0.5 15 3 10 1.E-08 1.E-12 19 500

Daniel J. Sharpe
April 2019
*/

#include "read_ktn.h"
#include "ktn.h"
#include "utils.h"
#include "quality.h"
#include <iostream>
#include <array>
#include <vector>
#include <algorithm>
#include <random>
#include <numeric>
#include <limits>
#include <map>
#include <iterator>
#include <stdlib.h>
#include <cstdlib>
#include <boost/python.hpp>
#include <Python.h>

using namespace std;

/* class for multi-level regularised Markov clustering (MLR-MCL) */
class MLR_MCL {

    typedef vector<map<int,int>> nodemap_vec;
    typedef pair<vector<pair<double,int>>,vector<int>> Csr_mtx; // matrix in CSR (or CSC) sparse format

    public:
    MLR_MCL(double,double,int,int,int,double,double,unsigned int,int,int);
    ~MLR_MCL();
    void run_mcl(Network&);
    static void calc_quality_metrics(Network&,int);

    double r; double b; double eps; double tau;
    int n_C; int n_cur; int max_it; unsigned int seed; int min_C; int min_comm_sz;

    private:
    void mcl_main_ops(Csr_mtx&,const Csr_mtx&);
    void coarsen_graph(Network&);
    void heavy_edge_matching(Network&,int);
    Csr_mtx get_matrix_exp(Network&);
    Csr_mtx get_init_sparse_mtx(Network&);
    Csr_mtx regularise(const Csr_mtx&,const Csr_mtx&);
    void inflate(Csr_mtx&);
    void prune(Csr_mtx&);
    void renormalise(Csr_mtx&);
    Csr_mtx get_reg_mtx(const Csr_mtx&);
    Csr_mtx project_flow(Network&,Csr_mtx&,const map<int,int>&);
    void interpret_clust(Network&,const Csr_mtx&);

    vector<int> g_C_size; // size of coarsened graph at each step
    nodemap_vec nodemap; // mapping of nodes to higher levels in multi-level coarsening procedure (pairwise merging)
    map<int,int> idxmap; // mapping of min ID's of coarsened graph to indices of corresponding transition matrix
    vector<int> idxlist; // ordered list of min ID's (posns in list correspond to indices of transition matrix)
    int n_row_occ; // number of rows of the stochastic transition matrix that are currently occupied
};

MLR_MCL::MLR_MCL(double d1,double d2,int i1,int i2,int i3,double d3,double d4, \
                 unsigned int i4,int i5,int i6) : \
                r(d1),b(d2),n_C(i1),n_cur(i2),max_it(i3),eps(d3),tau(d4),seed(i4), \
                min_C(i5), min_comm_sz(i6) {
    srand(seed);
    nodemap.resize(n_C);
};

MLR_MCL::~MLR_MCL() {};

/* main loop to drive multi-level regularised Markov clustering */
void MLR_MCL::run_mcl(Network &ktn) {
    coarsen_graph(ktn);
    Csr_mtx t_mtx_sp = get_matrix_exp(ktn); // transition matrix (CSR format)
    renormalise(t_mtx_sp);
    Csr_mtx tG_mtx_sp = get_reg_mtx(t_mtx_sp); // regularisation matrix (CSC format)
//    cout << "Initial matrix:" << endl; print_sparse_matrix(t_mtx_sp);
//    cout << "Initial reg matrix:" << endl; print_sparse_matrix(tG_mtx_sp);
    // run curtailed MLR-MCL
    for (int i=g_C_size.size()-1;i>=0;i--) {
        cout << ">>>>> running curtailed R-MCL on coarsened graph at level " << i+1 << \
                "  size of graph: " << t_mtx_sp.second.size() << endl;
        for (int j=0;j<n_cur;j++) { mcl_main_ops(t_mtx_sp,tG_mtx_sp); }
        cout << "calling project_flow()..." << endl;
        t_mtx_sp = project_flow(ktn,t_mtx_sp,nodemap[i]); // refined transition matrix
        cout << "returned projected t_mtx_sp. Renormalising..." << endl;
        renormalise(t_mtx_sp);
//        print_sparse_matrix(t_mtx_sp);
        tG_mtx_sp = get_reg_mtx(t_mtx_sp);
        cout << "returned regularisation matrix" << endl;
    }
    for (int i=0;i<max_it;i++) { cout << ">>>>> MCL iteration " << i << endl;
        mcl_main_ops(t_mtx_sp,tG_mtx_sp); }
    cout << ">>>>> End of MLR-MCL" << endl;
    interpret_clust(ktn,t_mtx_sp);
//    cout << "final matrix..." << endl; print_sparse_matrix(t_mtx_sp);
    if (min_comm_sz>0) Quality_clust::post_processing(ktn,min_comm_sz);
    cout << ">>>>> Writing communities and attractors to file..." << endl;
    Quality_clust::write_comms(ktn);
}

/* operations of Markov clustering main loop */
void MLR_MCL::mcl_main_ops(Csr_mtx &t_mtx_sp, const Csr_mtx &tG_mtx_sp) {
    t_mtx_sp = regularise(t_mtx_sp,tG_mtx_sp);
//    print_sparse_matrix(t_mtx_sp);
    inflate(t_mtx_sp);
//    print_sparse_matrix(t_mtx_sp);
    cout << "no. elems before prune " << t_mtx_sp.first.size() << " occ rows before prune " << n_row_occ << endl;
    prune(t_mtx_sp);
//    print_sparse_matrix(t_mtx_sp);
    renormalise(t_mtx_sp);
//    print_sparse_matrix(t_mtx_sp);
    cout << "no. elems after prune " << t_mtx_sp.first.size() << " occ rows after prune " << n_row_occ << endl;
}

/* calculate metrics to assess clustering quality. NB because HEM modifies the Network data
   structure, calling this function has to be done in a separate execution with the
   output_flag arg set to 1 */
void MLR_MCL::calc_quality_metrics(Network &ktn, int min_sz) {
    if (min_sz>0) Quality_clust::post_processing(ktn,min_sz);
    Quality_clust::find_intercomm_edges(ktn);
    cout << ">>>>> Writing inter-community edge bool values to file..." << endl;
    Quality_clust::find_intercomm_edges(ktn);
    cout << ">>>>> Modularity Q: " << Quality_clust::calc_modularity(ktn) << endl;
    double Ncut = Quality_clust::calc_avgncut(ktn);
    cout << ">>>>> Normalised cut: " << Ncut*double(ktn.n_comms) << endl;
    cout << ">>>>> Avg. Normalised cut: " << Ncut << endl;
    double conduc = Quality_clust::calc_conductance(ktn);
    cout << ">>>>> Conductance: " << conduc*double(ktn.n_comms) << endl;
    cout << ">>>>> Avg. conductance: " << conduc << endl;
}

/* heavy edge matching for graph coarsening */
void MLR_MCL::heavy_edge_matching(Network &ktn, int i_C) {
    Edge *edgeptr;
    vector<int> node_ids(ktn.tot_nodes);
    iota(begin(node_ids),end(node_ids),1);
    random_shuffle(begin(node_ids),end(node_ids));
    int matchnode_id;
    for (int i=0;i<ktn.tot_nodes;i++) {
        ktn.min_nodes[i].hem_flag = false; }
    for (int i=0;i<ktn.tot_nodes;i++) {
        // check if node has already been matched or deleted
        if ((ktn.min_nodes[node_ids[i]-1].hem_flag) || (ktn.min_nodes[node_ids[i]-1].deleted)) { continue; }
        else { ktn.min_nodes[node_ids[i]-1].hem_flag = true; }
        // find the neighbour TO node i associated with edge of greatest weight, this node will become subsumed
        edgeptr = ktn.min_nodes[node_ids[i]-1].top_to;
        double best_w = -numeric_limits<double>::infinity(); matchnode_id = 0;
        if (edgeptr != nullptr) {
        do {
            if ((edgeptr->from_node->hem_flag) || \
                (ktn.min_nodes[edgeptr->from_node->min_id-1].deleted)) {
                edgeptr = edgeptr->next_to; continue; }
            if (edgeptr->w > best_w) {
                best_w = edgeptr->w;
                matchnode_id = edgeptr->from_node->min_id;
            }
            edgeptr = edgeptr->next_to;
        } while (edgeptr!=nullptr);
        }
        else { cout << "Fatal error: node " << node_ids[i] << " is disconnected" << endl; exit(EXIT_FAILURE); }
        if (matchnode_id != 0) { // found a matching, merge and update nodemap
            ktn.min_nodes[matchnode_id-1].hem_flag = true;
            ktn.merge_nodes(ktn.min_nodes[node_ids[i]-1].min_id-1,matchnode_id-1);
            nodemap[i_C][node_ids[i]] = matchnode_id;
        } else { // node is not matched, maps to itself
            nodemap[i_C][node_ids[i]] = node_ids[i]; }
    }
}

/* function for coarsening of the graph */
void MLR_MCL::coarsen_graph(Network &ktn) {
    int i=0;
    if (ktn.n_nodes > min_C) {
    cout << "Coarsening the network..." << endl;
    do {
        heavy_edge_matching(ktn,i);
        i++;
        cout << ">>>> level " << i << ": n_nodes is now: " << ktn.n_nodes << endl;
        g_C_size.emplace_back(ktn.n_nodes);
    } while ((i < n_C) && (ktn.n_nodes > min_C));
    }
    if ((min_C > 0) && (ktn.n_nodes > min_C)) {
        cout << "Error: no. of nodes in most coarse graph exceeds maximum allowed" << endl; exit(EXIT_FAILURE); }
    n_row_occ = ktn.n_nodes;
    cout << "finished graph coarsening after " << i << " iterations" << endl;
}

/* call external Python script to compute the matrix exponential of the coarsened graph.
   Return the initial, coarsened column-stochastic flow matrix in CSR sparse format */
MLR_MCL::Csr_mtx MLR_MCL::get_matrix_exp(Network &ktn) {

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
    vector<int> k_rl(g_C_size.back(),0); // cumulative row lengths of transition rate matrix
    vector<double> k_elems; vector<int> k_minids; // temporary vectors for matrix elems and column indices
    Edge *edgeptr;
    // set the elements of the (flattened) transition rate matrix
    int rn=0;
    for (int i=0;i<ktn.tot_nodes;i++) {
        int rl=0;
        if (ktn.min_nodes[i].deleted) { continue; }
        edgeptr = ktn.min_nodes[i].top_to;
        bool edge_exist = false;
        if (edgeptr!=nullptr) {
        do {
            if (edgeptr->from_node->min_id-1==i) {
                cout << "Error: node " << i << " is connected to itself after coarsening" << endl;
                throw Network::Ktn_exception(); }
            // note that the column indices/min id's are not initially in order
            if (edgeptr->from_node->deleted) { edgeptr = edgeptr->next_to; continue; }
            if (!edge_exist) { edge_exist = true; }
            k_elems.emplace_back(edgeptr->w); k_minids.emplace_back(edgeptr->from_node->min_id);
            try {
                int dummy = idxmap.at(edgeptr->from_node->min_id);
            } catch (const out_of_range& oor) {
                idxmap[edgeptr->from_node->min_id] = -1; } // initial dummy value in map
            rl++;
            edgeptr = edgeptr->next_to;
        } while (edgeptr != nullptr);
        }
        if (!edge_exist) { cout << "Error: node " << i << " is disconnected after coarsening" << endl;
            throw Network::Ktn_exception(); }
        else { k_rl[rn]=rl; rn++; }
    }
    for (vector<int>::iterator it=k_rl.begin()+1;it!=k_rl.end();it++) {
        *it += *(it-1); }
    // initial setup of index list and index map
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

/* perform regularisation (matrix product) operation for the transition matrix */
MLR_MCL::Csr_mtx MLR_MCL::regularise(const Csr_mtx &T_csr, const Csr_mtx &TG_csr) {

    Csr_mtx T_csr_new;
    T_csr_new.first.reserve(n_row_occ*T_csr.second.size());
    T_csr_new.second.resize(T_csr.second.size());
    for (int i=0;i<T_csr.second.size();i++) { T_csr_new.second[i]=0; }
    signed int lo=0, hi=T_csr.second[0]; // low/high indices of elements of current row in transition matrix
    cout << "in regularise" << endl;
    int rn=0;
    int mop=0;
    const vector<pair<double,int>> *const vecptr1 = &T_csr.first, *const vecptr2 = &TG_csr.first;
    for (int i=0;i<n_row_occ;i++) {
        while (lo==hi) {
            lo = hi; hi = T_csr.second[rn+1]; rn++; }
        // low/high indices of elements of current col in regularisation matrix
        signed int lo_col=0, hi_col=TG_csr.second[0];
        int cn=0;
        for (int j=0;j<T_csr.second.size();j++) {
            if (lo_col==hi_col) {
                cout << "Error: column " << j+1 << " of regularisation matrix is empty" << endl;
                exit(EXIT_FAILURE); }
            int curr_row=lo, curr_col=lo_col;
            double val=0.;
            // multiplication of two sparse vectors with O(m+n) complexity
            while ((curr_row < hi) && (curr_col < hi_col)) {
                if ((*vecptr1)[curr_row].second < (*vecptr2)[curr_col].second) { curr_row++;
                } else if ((*vecptr1)[curr_row].second > (*vecptr2)[curr_col].second) { curr_col++;
                } else {
                    val += (*vecptr1)[curr_row].first*(*vecptr2)[curr_col].first;
                    curr_row++; curr_col++; mop++;
                }
            }
            if (val>eps) { T_csr_new.first.emplace_back(make_pair(val,cn)); T_csr_new.second[rn]++; }
            lo_col = hi_col; hi_col = TG_csr.second[cn+1]; cn++;
        }
        lo = hi; hi = T_csr.second[rn+1]; rn++;
    }
    for (int i=1;i<T_csr_new.second.size();i++) {
        T_csr_new.second[i] += T_csr_new.second[i-1]; }
    T_csr_new.first.shrink_to_fit();
    cout << "number of mop's: " << mop << endl;
    cout << "leaving regularise()" << endl;
    return T_csr_new;
}

/* inflation operation for the transition matrix */
void MLR_MCL::inflate(Csr_mtx &T_csr) {
    for (int i=0;i<T_csr.first.size();i++) {
        T_csr.first[i].first = pow(T_csr.first[i].first,r); }
}

/* prune small values from sparse transition matrix */
void MLR_MCL::prune(Csr_mtx &T_csr) {

    vector<int> T_rl = T_csr.second; // raw row lengths
    for (vector<int>::iterator it=T_rl.end()-1;it!=T_rl.begin();it--) {
        *it -= *(it-1); }
    vector<pair<double,int>>::iterator it_vec, new_end;
    int k=0, rn=0; if (k==T_csr.second[rn]) while (k==T_csr.second[rn]) { rn++; };
    new_end = remove_if(T_csr.first.begin(),T_csr.first.end(),[=,&k,&rn,&T_rl](pair<double,int> ec) {
        if (ec.first<eps) T_rl[rn]--;
        k++;
        if ((k==T_csr.second[rn]) && (k<T_csr.first.size())) while (k==T_csr.second[rn]) { rn++; };
        return (ec.first<eps); });
    T_csr.first.erase(new_end,T_csr.first.end());
    int row_occ_count;
    T_csr.second[0] = T_rl[0];
    if (T_rl[0]==0) { row_occ_count=0; } else { row_occ_count=1; }
    for (int i=1;i<T_rl.size();i++) { // set new row lengths of transition matrix
        T_csr.second[i] = T_rl[i] + T_csr.second[i-1];
        if (T_rl[i]!=0) row_occ_count++;
    }
    n_row_occ = row_occ_count;
    if (T_csr.second.back()!=T_csr.first.size()) {
        cout << "Fatal error: row length vector of sparse matrix is incorrect" << endl; exit(EXIT_FAILURE); }
}

/* renormalise columns of sparse transition matrix */
void MLR_MCL::renormalise(Csr_mtx &T_csr) {

    vector<double> cum_sum(T_csr.second.size(),0.);
    vector<pair<double,int>>::iterator it_vec;
    for (it_vec=T_csr.first.begin();it_vec!=T_csr.first.end();it_vec++) { // accumulate column sums
        cum_sum[it_vec->second] += it_vec->first;
    }
    for (it_vec=T_csr.first.begin();it_vec!=T_csr.first.end();it_vec++) { // renormalise columns
        it_vec->first *= (1./cum_sum[it_vec->second]); }
}

/* given a transition matrix in CSR format, return a regularisation matrix in CSC format */
MLR_MCL::Csr_mtx MLR_MCL::get_reg_mtx(const Csr_mtx &T_csr) {

    cout << "calculating regularisation matrix..." << endl;
    vector<vector<pair<double,int>>> TG_er_al(T_csr.second.size()); // elems and row indices, adj list
    vector<int> TG_cl(T_csr.second.size()); // col lengths
    vector<pair<double,int>>::const_iterator it_vec;
    int k=0, rn=0; if (k==T_csr.second[rn]) while (k==T_csr.second[rn]) { rn++; };
    for (it_vec=T_csr.first.begin();it_vec!=T_csr.first.end();it_vec++) {
        TG_er_al[it_vec->second].emplace_back(make_pair(it_vec->first,rn));
        k++;
        if ((k==T_csr.second[rn]) && (k<T_csr.first.size())) while (k==T_csr.second[rn]) { rn++; };
    }
    k=0; TG_cl[k] = TG_er_al[k].size();
    for (const auto &colvec: TG_er_al) {
        if (colvec.size()==0) {
            cout << "Error: column " << k+1 << " is unoccupied" << endl; exit(EXIT_FAILURE); }
        if (k>=1) TG_cl[k] = TG_cl[k-1] + colvec.size();
        k++; }
    vector<pair<double,int>> TG_er = flatten<pair<double,int>>(TG_er_al);
    Csr_mtx TG_csr = make_pair(TG_er,TG_cl);
    return TG_csr;
}

/* use nodemaps to undo the coarsening of the graph (refinement) */
MLR_MCL::Csr_mtx MLR_MCL::project_flow(Network &ktn, Csr_mtx &T_csr, const map<int,int> &curr_nodemap) {

    vector<int> idxlist_old = idxlist;
    map<int,int>::const_iterator it_map;
    for (it_map=curr_nodemap.begin();it_map!=curr_nodemap.end();it_map++) { // update the list of indices
        if (it_map->second != it_map->first) { // add second node of pair to list of indices
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
    int rn_occ=0; // occupied row no.
    if (k==T_csr.second[rn]) while (k==T_csr.second[rn]) { rn++; }; // account for if first elem is not in first row
    vector<pair<double,int>>::iterator it_vec;
    // only two of the four corresponding elements in the refined matrix are non-zero
    for (it_vec=T_csr.first.begin();it_vec!=T_csr.first.end();it_vec++) {
        T_ec_new[idxmap.at(idxlist_old[rn])].emplace_back(make_pair(it_vec->first, \
                 idxmap.at(idxlist_old[it_vec->second])));
        if (idxlist_old[it_vec->second] != curr_nodemap.at(idxlist_old[it_vec->second])) {
        T_ec_new[idxmap[idxlist_old[rn]]].emplace_back(make_pair(it_vec->first, \
                 idxmap.at(curr_nodemap.at(idxlist_old[it_vec->second]))));
        }
        k++;
        if ((k==T_csr.second[rn]) && (k<T_csr.first.size())) while (k==T_csr.second[rn]) { rn++; };
    }
    k=0;
    for (auto &rowvec: T_ec_new) {
        sort(rowvec.begin(),rowvec.end(),[](pair<double,int> pair1,pair<double,int> pair2) {
            return (pair1.second < pair2.second); });
        T_rl[k] = rowvec.size(); if (k>0) T_rl[k] += T_rl[k-1];
        k++; }
    // flatten the refined transition matrix
    vector<pair<double,int>> T_ec = flatten<pair<double,int>>(T_ec_new);
    cout << "no. of elems of matrix is now: " << T_ec.size() << endl;
    Csr_mtx T_csr_new = make_pair(T_ec,T_rl);
    return T_csr_new;
}

/* after the clustering procedure has finished, interpret the flow matrix as a clustering
   characterised by attractors */
void MLR_MCL::interpret_clust(Network &ktn, const Csr_mtx &T_csr) {

    cout << ">>>>> Interpreting the final stochastic matrix as a clustering..." << endl;
    int n_comm=-1; // counter for community IDs
    vector<pair<double,int>>::const_iterator it_vec;
    int k=0, rn=0; if (k==T_csr.second[rn]) while (k==T_csr.second[rn]) { rn++; };
    for (it_vec=T_csr.first.begin();it_vec!=T_csr.first.end();it_vec++) {
        if (it_vec->first > 0.1) { // value considered non-negligible
        if (!ktn.min_nodes[rn].attractor) {
            ktn.min_nodes[rn].attractor=true;
            if (ktn.min_nodes[rn].comm_id==-1) { n_comm++;
            } else {
                cout << "Warning: community " << ktn.min_nodes[rn].comm_id << \
                        " is characterised by more than one attractor" << endl; }
        }
        if (ktn.min_nodes[it_vec->second].comm_id==-1) ktn.min_nodes[it_vec->second].comm_id = n_comm;
        }
        k++;
        if (k==T_csr.second[rn]) while (k==T_csr.second[rn]) { rn++; };
    }
    ktn.n_comms = n_comm+1;
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
    unsigned int seed = stoi(argv[10]); // random seed
    int min_C; // min. no. of nodes in coarsened graph (optional)
    int min_comm_sz; // min. no. of (do post-processing if >0) (optional)
    int debug_flag; // run debug tests Y/N (optional)
    int output_flag; // read communities from file and calculate quality metrics only (optional)
    if (argc > 11) { min_C = stoi(argv[11]); } else { min_C = 0; }
    if (argc > 12) { min_comm_sz = stoi(argv[12]); } else { min_comm_sz = 0; }
    if (argc > 13) { debug_flag = stoi(argv[13]); } else { debug_flag = 0; }
    if (argc > 14) { output_flag = stoi(argv[14]); } else { output_flag = 0; }
    cout << ">>>>> Finished reading input arguments" << endl;

    vector<pair<int,int>> ts_conns = Read_ktn::read_double_col<int>(nts,"ts_conns.dat");
    vector<double> ts_weights = Read_ktn::read_single_col<double>(2*nts,"ts_weights.dat");
    vector<double> stat_probs = Read_ktn::read_single_col<double>(nmin,"stat_prob.dat");

    Network ktn(nmin,nts);
    Network::setup_network(ktn,nmin,nts,ts_conns,ts_weights,stat_probs);
    ts_conns.resize(0); ts_weights.resize(0); stat_probs.resize(0);

    if (debug_flag) { run_debug_tests(ktn); exit(0); }
    if (output_flag) {
        Quality_clust::read_comms(ktn);
        MLR_MCL::calc_quality_metrics(ktn,min_comm_sz); exit(0); }

    MLR_MCL mcl_obj (r,b,n_C,n_cur,max_it,eps,tau,seed,min_C,min_comm_sz);
    mcl_obj.run_mcl(ktn);

    cout << ">>>>> Finished" << endl;

    return 0;
}
