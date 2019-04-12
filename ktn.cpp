#include "ktn.h"
#include <vector>
#include <algorithm>

#include <iostream>

using namespace std;

Network::Network(int nmin, int nts) {
    min_nodes.resize(nmin);
    ts_edges.resize(2*nts);
    n_nodes = 0; n_edges = 0;
}

Network::~Network() {}

// delete node i 
void Network::del_node(int i) {
    if (min_nodes[i].deleted) { throw Ktn_exception(); }
    Edge *edgeptr;
    edgeptr = min_nodes[i].top_to;
    do {
        del_to_edge(i);
        edgeptr = edgeptr->next_to;
    } while (edgeptr != nullptr);
    edgeptr = min_nodes[i].top_from;
    do {
        del_from_edge(i);
        edgeptr = edgeptr->next_from;
    } while (edgeptr != nullptr);
    min_nodes[i].deleted = true;
    n_nodes--;
}

// edge j goes TO node i
void Network::add_to_edge(int i, int j) {
    if (min_nodes[i].top_to != nullptr) {
        ts_edges[j].next_to = min_nodes[i].top_to; }
    min_nodes[i].top_to = &ts_edges[j];
    n_edges++;
}

// edge j goes FROM node i
void Network::add_from_edge(int i, int j) {
    if (min_nodes[i].top_from != nullptr) {
        ts_edges[j].next_from = min_nodes[i].top_from; }
    min_nodes[i].top_from = &ts_edges[j];
    n_edges++;
}

// delete the top TO edge for node i
void Network::del_to_edge(int i) {
    if (min_nodes[i].top_to != nullptr) {
        if (min_nodes[i].top_to->next_to != nullptr) {
            min_nodes[i].top_to = min_nodes[i].top_to->next_to;
        } else {
            min_nodes[i].top_to = nullptr;
        }
        n_edges--;
    } else {
        throw Ktn_exception();
    }
}

// delete the top FROM edge for node i
void Network::del_from_edge(int i) {
    if (min_nodes[i].top_from != nullptr) {
        if (min_nodes[i].top_from->next_from != nullptr) {
            min_nodes[i].top_from = min_nodes[i].top_from->next_from;
        } else {
            min_nodes[i].top_from = nullptr;
        }
        n_edges--;
    } else {
        throw Ktn_exception();
    }
}

// merge a pair of nodes. The new sets of TO and FROM edges are the unions of all the TO and FROM edges
// for the constituent pair of nodes. cf node contraction operation. Node j is merged into node i.
void Network::merge_nodes(int i, int j) {
    if ((min_nodes[i].deleted) || (min_nodes[j].deleted)) { throw Ktn_exception(); }
    Edge *edgeptr;
    // get all TO and FROM neighbours for node i (need to check these against TO and FROM neighbours for node j)
    vector<int> node_i_nbrs_to; vector<int> node_i_nbrs_from;
    vector<int> ts_nbrs_to; vector<int> ts_nbrs_from;
    edgeptr = min_nodes[i].top_to;
    cout << "   is null? " << (edgeptr==nullptr) << endl;
    cout << "   to nodes..." << endl;
    do {
        cout << "      " << edgeptr->ts_id << endl;
        node_i_nbrs_to.push_back(edgeptr->to_node->min_id);
        ts_nbrs_to.push_back(edgeptr->ts_pos);
        edgeptr = edgeptr->next_to;
    } while (edgeptr != nullptr);
    edgeptr = min_nodes[i].top_from;
    cout << "   from nodes..." << endl;
    do {
        node_i_nbrs_from.push_back(edgeptr->from_node->min_id);
        ts_nbrs_from.push_back(edgeptr->ts_pos);
        edgeptr = edgeptr->next_from;
    } while (edgeptr != nullptr);
    // merge TO edges of j to node i
    cout << "   now do the merge" << endl;
    edgeptr = min_nodes[j].top_to;
    do {
        if (find(node_i_nbrs_to.begin(),node_i_nbrs_to.end(),edgeptr->to_node->min_id) != \
            node_i_nbrs_to.end()) {
// need to know posn of corresponding ts_pos entry in ts_nbrs_to
//            ts_edges[edgeptr->ts_id].w += ;
        } else {
            add_to_edge(i,edgeptr->ts_pos); }
        edgeptr = edgeptr->next_to;
    } while (edgeptr != nullptr);
    // merge FROM edges of j to node i
    // ...
    del_node(j);
}

// set up the kinetic transition network
void Network::setup_network(Network& ktn, int nmin, int nts, vector<pair<int,int>> ts_conns, \
                   vector<double> ts_weights) {

    ktn.n_nodes = nmin;
    for (int i=0;i<nmin;i++) {
        ktn.min_nodes[i].min_id = i+1;
    }
    for (int i=0;i<nts;i++) {
        ktn.ts_edges[2*i].ts_id = i+1;
        ktn.ts_edges[(2*i)+1].ts_id = i+1;
        ktn.ts_edges[2*i].ts_pos = 2*i;
        ktn.ts_edges[(2*i)+1].ts_id = (2*i)+1;
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

