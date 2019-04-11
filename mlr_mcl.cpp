/*
C++ code for multi-level regularised Markov clustering (MLR-MCL) of a kinetic transition network
Need two files: "ts_conns.dat" (two-column format: min1, min2) and "ts_weights.dat" (one-column format)
in the current directory

Compile with:
g++ -std=c++11 mlr_mcl.cpp read_ktn.cpp -o mlr_mcl

Execute with e.g.:
./mlr_mcl 22659 34145 1.5 0.5 10 4 1.D-06 1.D-06

Daniel J. Sharpe
April 2019
*/


#include "read_ktn.h"
#include <iostream>
#include <string>
#include <array>
#include <vector>

using namespace std;

struct Node;

struct Edge {
    int ts_id;
    double w;
    bool deadts = false;
    Node *to_node;
    Node *from_node;
    Edge *next_to;
    Edge *next_from;
};

struct Node {
    int min_id;
    Edge *top_to;
    Edge *top_from;
};

// structure containing the kinetic transition network
struct Network {
    Network(int,int);
    void add_to_edge(int,int);
    void add_from_edge(int,int);
    vector<Node> min_nodes;
    vector<Edge> ts_edges;
};

Network::Network(int nmin, int nts) {
    min_nodes.resize(nmin);
    ts_edges.resize(2*nts);
}

// edge j goes TO node i
void Network::add_to_edge(int i, int j) {
    if (min_nodes[i].top_to != nullptr) {
        ts_edges[j].next_to = min_nodes[i].top_to; }
    min_nodes[i].top_to = &ts_edges[j];
}

// edge j goes FROM node i
void Network::add_from_edge(int i, int j) {
    if (min_nodes[i].top_from != nullptr) {
        ts_edges[j].next_from = min_nodes[i].top_from; }
    min_nodes[i].top_from = &ts_edges[j];
}

// set up the kinetic transition network
void setup_network(Network& ktn, int nmin, int nts, vector<pair<int,int>> ts_conns, \
                   vector<double> ts_weights) {

    for (int i=0;i<nmin;i++) {
        ktn.min_nodes[i].min_id = i+1;
    }
    for (int i=0;i<nts;i++) {
        ktn.ts_edges[2*i].ts_id = i+1;
        ktn.ts_edges[(2*i)+1].ts_id = i+1;
        if (ts_conns[i].first == ts_conns[i].second) {
            ktn.ts_edges[2*i].deadts = true;
            ktn.ts_edges[(2*i)+1].deadts = true;
            continue; }
        ktn.ts_edges[2*i].w = ts_weights[2*i];
        ktn.ts_edges[(2*i)+1].w = ts_weights[(2*i)+1];
        ktn.ts_edges[2*i].from_node = &ktn.min_nodes[ts_conns[i].first-1];
        ktn.ts_edges[2*i].to_node = &ktn.min_nodes[ts_conns[i].second-1];
        ktn.ts_edges[(2*i)+1].from_node = &ktn.min_nodes[ts_conns[i].second-1];
        ktn.ts_edges[(2*i)+1].to_node = &ktn.min_nodes[ts_conns[i].first-1];

        ktn.add_to_edge(ts_conns[i].second-1,2*i);
        ktn.add_from_edge(ts_conns[i].first-1,2*i);
        ktn.add_to_edge(ts_conns[i].first-1,(2*i)+1);
        ktn.add_from_edge(ts_conns[i].second-1,(2*i)+1);
    }
}

// class for multi-level regularised Markov clustering (MLR-MCL)
class MLR_MCL {

    public:
    MLR_MCL(double,double,int,int,double,double);
    ~MLR_MCL();
    void run_mcl();

    double r; double b; double eps; double tau;
    int n_C; int n_cur;

    private:
    void coarsen_graph();
    void construct_matrix();
    void prune();
    void regularise();
    void expand();
    void inflate();
};

MLR_MCL::MLR_MCL(double d1,double d2,int i1,int i2,double d3, double d4) : \
                r(d1),b(d2),n_C(i1),n_cur(i2),eps(d3),tau(d4) {};
MLR_MCL::~MLR_MCL() {};

void MLR_MCL::run_mcl() {

}

// heavy edge matching for graph coarsening
//pair<vector<int>,vector<int>> heavy_edge_matching() {
//}

// recursive function for coarsening of the graph
void MLR_MCL::coarsen_graph() {

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

    vector<pair<int,int>> ts_conns = Read_ktn::read_ts_conns(nts);
    vector<double> ts_weights = Read_ktn::read_ts_weights(nts);

    Network ktn(nmin,nts);
    setup_network(ktn,nmin,nts,ts_conns,ts_weights);

    // TESTS
    cout << ktn.min_nodes[3].min_id << " (should be 4)" << endl;
    cout << ktn.ts_edges[2*500].from_node->min_id << " " << ktn.ts_edges[2*500].to_node->min_id << endl;
    cout << ktn.ts_edges[(2*500)+1].from_node->min_id << " " << ktn.ts_edges[(2*500)+1].to_node->min_id << endl;
    cout << "iterate over all edges pointing TO node 1..." << endl;
    int i=0;
    Node *nodeptr; Edge *edgeptr;
    edgeptr = ktn.min_nodes[i].top_to;
    do {
        cout << edgeptr->ts_id << "  " << edgeptr->w << "    " << edgeptr->from_node->min_id << " " << edgeptr->to_node->min_id << endl;
        edgeptr = edgeptr->next_to;
    } while (edgeptr != nullptr);
    cout << "iterate over all edges pointing FROM node 1..." << endl;
    edgeptr = ktn.min_nodes[i].top_from;
    do {
        cout << edgeptr->ts_id << "  " << edgeptr->w << "    " << edgeptr->from_node->min_id << " " << edgeptr->to_node->min_id << endl;
        edgeptr = edgeptr->next_from;
    } while (edgeptr != nullptr);

    ts_conns.resize(0); ts_weights.resize(0);


    MLR_MCL mcl_obj (r,b,n_C,n_cur,eps,tau);
}
