/*
C++ code for multi-level regularised Markov clustering (MLR-MCL) of a kinetic transition network
Need two files: "ts_conns.dat" (two-column format: min1, min2) and "ts_weights.dat" (one-column format)
in the current directory

Compile with:
g++ -std=c++11 mlr_mcl.cpp ktn.cpp read_ktn.cpp utils.cpp -o mlr_mcl

Execute with e.g.:
./mlr_mcl 22659 34145 1.5 0.5 10 4 1.D-06 1.D-06 19

Daniel J. Sharpe
April 2019
*/

#include "read_ktn.h"
#include "ktn.h"
#include "utils.h"
#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>

using namespace std;

// class for multi-level regularised Markov clustering (MLR-MCL)
class MLR_MCL {

    typedef vector<vector<pair<int,int>>> nodemap;

    public:
    MLR_MCL(double,double,int,int,double,double,int,int);
    ~MLR_MCL();
    void run_mcl(Network&);

    double r; double b; double eps; double tau;
    int n_C; int n_cur; int seed; int min_C;

    private:
    void coarsen_graph(Network&);
    void heavy_edge_matching(Network&);
    void construct_matrix();
    void prune();
    void regularise();
    void expand();
    void inflate();

    nodemap nodemap1; nodemap nodemap2;
};

MLR_MCL::MLR_MCL(double d1,double d2,int i1,int i2,double d3, double d4, int i3, int i4) : \
                r(d1),b(d2),n_C(i1),n_cur(i2),eps(d3),tau(d4),seed(i3),min_C(i4) {

    nodemap1.resize(n_C); nodemap2.resize(n_C);
};

MLR_MCL::~MLR_MCL() {};

void MLR_MCL::run_mcl(Network &ktn) {
    coarsen_graph(ktn);
}

// heavy edge matching for graph coarsening
void MLR_MCL::heavy_edge_matching(Network &ktn) {
    Edge *edgeptr;
    vector<int> node_ids(ktn.tot_nodes);
    iota(begin(node_ids),end(node_ids),1);
    random_shuffle(begin(node_ids),end(node_ids));
    int matchnode_id;
    for (int i=0;i<ktn.tot_nodes;i++) {
        ktn.min_nodes[i].hem_flag = false; }
    for (int i=0;i<ktn.tot_nodes;i++) {
        cout << " iteration no. " << i+1 << " chosen node... " << ktn.min_nodes[node_ids[i]-1].min_id << endl;
        // check if node has already been matched or deleted
        if ((ktn.min_nodes[node_ids[i]-1].hem_flag) || (ktn.min_nodes[node_ids[i]-1].deleted)) { continue; }
        else { ktn.min_nodes[node_ids[i]-1].hem_flag = true; }
//        cout << "  has not been deleted or matched..." << endl;
        // find the neighbour FROM node i with largest weight
        edgeptr = ktn.min_nodes[node_ids[i]-1].top_from;
//        cout << "edgeptr FROM " << edgeptr->from_node->min_id << " edgeptr TO " << edgeptr->to_node->min_id << endl;
        double best_w = -numeric_limits<double>::infinity(); matchnode_id = 0;
        if (edgeptr != nullptr) {
        do {
//            cout << "nbr node " << edgeptr->to_node->min_id << " w " << edgeptr->w << endl;
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
/*
        if (node_ids[i] == 21625) {
            cout << "    TO nodes for node " << node_ids[i] << "..." << endl;
            Edge *edgeptr = ktn.min_nodes[node_ids[i]-1].top_to;
            if (edgeptr!=nullptr) { print_edgeptr_info(edgeptr,1);
            } else { cout << "    no top TO node! " << endl; }
            cout << "    FROM nodes for node " << node_ids[i] << "..." << endl;
            edgeptr = ktn.min_nodes[node_ids[i]-1].top_from;
            if (edgeptr!=nullptr) { print_edgeptr_info(edgeptr,2);
            } else { cout << "    no top FROM node! " << endl; }
        }
*/
        if (matchnode_id != 0) { // found a matching
            cout << "    matched node been deleted? " << ktn.min_nodes[matchnode_id-1].deleted << \
                    " or matched? " << ktn.min_nodes[matchnode_id-1].hem_flag << endl;
            cout << "    highest weight FROM nbr is node " << matchnode_id << " w " << best_w << endl;
            ktn.min_nodes[matchnode_id-1].hem_flag = true;
            ktn.merge_nodes(node_ids[i]-1,matchnode_id-1);
            // update nodemaps
        } else { ktn.min_nodes[node_ids[i]-1].deleted = true; ktn.n_nodes--; }
//        break;
    }
}

// function for coarsening of the graph
void MLR_MCL::coarsen_graph(Network &ktn) {
    int i=0;
    do {
        heavy_edge_matching(ktn);
        i++;
        cout << ">>>> n_nodes IS NOW: " << ktn.n_nodes << endl;
        break;
    } while ((i < n_C) && (ktn.n_nodes > min_C));
}

// construct a sparse matrix for matrix multiplication
void MLR_MCL::construct_matrix() {

}

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
