#include "ktn.h"
#include <vector>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <iostream>

using namespace std;

Network::Network(int nmin, int nts) {
    min_nodes.resize(nmin);
    ts_edges.resize(2*nts);
    n_nodes = 0; n_edges = 0; tot_nodes = 0; tot_edges = 0;
}

Network::~Network() {}

// delete node i 
void Network::del_node(int i) {
    if (min_nodes[i].deleted) { throw Ktn_exception(); }
    Edge *edgeptr;
    edgeptr = min_nodes[i].top_to;
    if (edgeptr!=nullptr) {
    do {
        del_to_edge(i);
        edgeptr = edgeptr->next_to;
    } while (edgeptr != nullptr);
    }
    edgeptr = min_nodes[i].top_from;
    if (edgeptr!=nullptr) {
    do {
        del_from_edge(i);
        edgeptr = edgeptr->next_from;
    } while (edgeptr != nullptr);
    }
    min_nodes[i].deleted = true;
    n_nodes--;
}

// edge j goes TO node i
void Network::add_to_edge(int i, int j) {
    if (min_nodes[i].top_to != nullptr) {
        ts_edges[j].next_to = min_nodes[i].top_to;
        min_nodes[i].top_to = &ts_edges[j]; }
    else {
        min_nodes[i].top_to = &ts_edges[j];
        min_nodes[i].top_to->next_to = nullptr; }
    n_edges++;
}

// edge j goes FROM node i
void Network::add_from_edge(int i, int j) {
    if (min_nodes[i].top_from != nullptr) {
        ts_edges[j].next_from = min_nodes[i].top_from;
        min_nodes[i].top_from = &ts_edges[j]; }
    else {
        min_nodes[i].top_from = &ts_edges[j];
        min_nodes[i].top_from->next_from = nullptr; }
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

// delete TO edge with ts_id j for node i
void Network::del_spec_to_edge(int i, int j) {
    Edge *edgeptr; Edge *edgeptr_prev = nullptr;
    bool ts_exists = false;
    if (min_nodes[i].top_to != nullptr) {
        edgeptr = min_nodes[i].top_to;
        do {
            if (edgeptr->ts_id==j) {
                if (edgeptr_prev != nullptr) {
                    edgeptr_prev->next_to = edgeptr->next_to;
                } else if (edgeptr->next_to != nullptr) {
                    min_nodes[i].top_to = edgeptr->next_to;
                } else if (edgeptr->next_to == nullptr) {
                    min_nodes[i].top_to = nullptr;
                }
                ts_exists = true; break;
            }
            edgeptr_prev = edgeptr;
            edgeptr = edgeptr->next_to;
        } while (edgeptr != nullptr);
    } else {
        cout << "Error: no TO edge to node with min_id " << i+1 << endl;
        throw Ktn_exception();
    }
    if (!ts_exists) { cout << "ts does not exist!" << endl; throw Ktn_exception(); }
}

// delete FROM edge with ts_id j for node i
void Network::del_spec_from_edge(int i, int j) {
    Edge *edgeptr; Edge *edgeptr_prev = nullptr;
    bool ts_exists = false;
    if (min_nodes[i].top_from != nullptr) {
        edgeptr = min_nodes[i].top_from;
        do {
            if (edgeptr->ts_id==j) {
                if (edgeptr_prev != nullptr) {
                    edgeptr_prev->next_from = edgeptr->next_from;
                } else if (edgeptr->next_from != nullptr) {
                    min_nodes[i].top_from = edgeptr->next_from;
                } else if (edgeptr->next_from ==nullptr) {
                    min_nodes[i].top_from = nullptr;
                }
                ts_exists = true; break;
            }
            edgeptr_prev = edgeptr;
            edgeptr = edgeptr->next_from;
        } while (edgeptr != nullptr);
    } else {
        cout << "Error: no FROM edge to node with min_id " << i+1 << endl;
        throw Ktn_exception();
    }
    if (!ts_exists) { throw Ktn_exception(); }
}

// update edge so that it now points TO i
void Network::update_to_edge(int i, int j) {
    Edge *edgeptr;
    edgeptr = &ts_edges[j];
    int old_to = edgeptr->to_node->min_id;
    edgeptr->to_node = &min_nodes[i];
    del_spec_to_edge(old_to-1,edgeptr->ts_id);
    add_to_edge(i,j);
}

// update edge so that it now points FROM i
void Network::update_from_edge(int i, int j) {
    Edge *edgeptr;
    edgeptr = &ts_edges[j];
    int old_from = edgeptr->from_node->min_id;
    edgeptr->from_node = &min_nodes[i];
    del_spec_from_edge(old_from-1,edgeptr->ts_id);
    add_from_edge(i,j);
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
    if (edgeptr != nullptr) {
    do {
        if (edgeptr->from_node->min_id==j+1) { del_spec_to_edge(i,edgeptr->ts_id);
            edgeptr = edgeptr->next_to; continue; }
        node_i_nbrs_to.push_back(edgeptr->from_node->min_id);
        ts_nbrs_to.push_back(edgeptr->ts_pos);
        edgeptr = edgeptr->next_to;
    } while (edgeptr != nullptr);
    }
    edgeptr = min_nodes[i].top_from;
    if (edgeptr != nullptr) {
    do {
        if (edgeptr->to_node->min_id==j+1) { del_spec_from_edge(i,edgeptr->ts_id);
            edgeptr = edgeptr->next_from; continue; }
        node_i_nbrs_from.push_back(edgeptr->to_node->min_id);
        ts_nbrs_from.push_back(edgeptr->ts_pos);
        edgeptr = edgeptr->next_from;
    } while (edgeptr != nullptr);
    }
    // merge TO edges of j to node i
    vector<int> to_edges_toadd; vector<int> from_edges_toadd;
    edgeptr = min_nodes[j].top_to;
    if (edgeptr != nullptr) {
    do {
        if (edgeptr->from_node->min_id==i+1 || edgeptr->from_node->deleted) { edgeptr = edgeptr->next_to; continue; }
        vector<int>::iterator it_find = find(node_i_nbrs_to.begin(),node_i_nbrs_to.end(), \
            edgeptr->from_node->min_id);
        if (it_find != node_i_nbrs_to.end()) { // this node TO j is also already a node TO i
            int nbr_idx = distance(node_i_nbrs_to.begin(),it_find);
            ts_edges[ts_nbrs_to[nbr_idx]].w = log(exp(ts_edges[ts_nbrs_to[nbr_idx]].w) + exp(ts_edges[edgeptr->ts_pos].w));
            del_spec_to_edge(edgeptr->from_node->min_id-1,edgeptr->ts_id);
        } else { // this node TO j does not already exist TO i
            to_edges_toadd.emplace_back(edgeptr->ts_pos);
        }
        edgeptr = edgeptr->next_to;
    } while (edgeptr != nullptr);
    }
    for (auto ts_pos: to_edges_toadd) {
        update_to_edge(i,ts_pos); }
    // merge FROM edges of j to node i
    edgeptr = min_nodes[j].top_from;
    if (edgeptr != nullptr) {
    do {
        if (edgeptr->to_node->min_id==i+1 || edgeptr->to_node->deleted) { edgeptr = edgeptr->next_from; continue; }
        vector<int>::iterator it_find = find(node_i_nbrs_from.begin(),node_i_nbrs_from.end(), \
            edgeptr->to_node->min_id);
        if (it_find != node_i_nbrs_from.end()) {
            int nbr_idx = distance(node_i_nbrs_from.begin(),it_find);
            ts_edges[ts_nbrs_from[nbr_idx]].w = log(exp(ts_edges[ts_nbrs_from[nbr_idx]].w) + exp(ts_edges[edgeptr->ts_pos].w));
            del_spec_from_edge(edgeptr->to_node->min_id-1,edgeptr->ts_id);
        } else {
            from_edges_toadd.emplace_back(edgeptr->ts_pos);
        }
        edgeptr = edgeptr->next_from;
    } while (edgeptr != nullptr);
    }
    for (auto ts_pos: from_edges_toadd) {
        update_from_edge(i,ts_pos); }
    // now delete node j from the network
    del_node(j);
}

/* calculate the in-degree or out_degree for node i */
double Network::calc_deg_inout(const Network &ktn, int i, int inout) {
    Edge *edgeptr;
    double deg_inout=0.;
    if (inout==1) { edgeptr = ktn.min_nodes[i].top_to; }
    else if (inout==2) { edgeptr = ktn.min_nodes[i].top_from; }
    if (edgeptr!=nullptr) {
    do {
        deg_inout = log(exp(deg_inout) + exp(edgeptr->w));
        if (inout==1) { edgeptr = edgeptr->next_to; }
        else if (inout==2) { edgeptr = edgeptr->next_from; }
    } while (edgeptr!=nullptr);
    }
    return deg_inout;
}

/* set up the kinetic transition network */
void Network::setup_network(Network& ktn, int nmin, int nts, const vector<pair<int,int>> ts_conns, \
                   const vector<double> ts_weights, const vector<double> stat_probs) {

    ktn.n_nodes = nmin; ktn.tot_nodes = nmin; ktn.n_edges = 2*nts; ktn.tot_edges = 2*nts;
    double tot_peq = 0.;
    for (int i=0;i<nmin;i++) {
        ktn.min_nodes[i].min_id = i+1;
        ktn.min_nodes[i].peq = stat_probs[i];
        tot_peq = log(exp(tot_peq) + exp(stat_probs[i]));
    }
    if (abs(tot_peq-1.)>1.E-15) {
        cout << "Error: total equilibrium probabilities of minima is: " << tot_peq << " =/= 1" << endl;
        throw Network::Ktn_exception(); }
    for (int i=0;i<nts;i++) {
        ktn.ts_edges[2*i].ts_id = i+1;
        ktn.ts_edges[(2*i)+1].ts_id = i+1;
        ktn.ts_edges[2*i].ts_pos = 2*i;
        ktn.ts_edges[(2*i)+1].ts_pos = (2*i)+1;
        if (ts_conns[i].first == ts_conns[i].second) {
            ktn.ts_edges[2*i].deadts = true;
            ktn.ts_edges[(2*i)+1].deadts = true;
            ktn.n_dead++;
            ktn.n_edges = ktn.n_edges-2;
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
    cout << "Finished reading in kinetic transition network" << endl;
    // account for if there are 2 edges connecting the same pair of nodes
    Edge *edgeptr;
    for (int i=0;i<nmin;i++) {
        edgeptr = ktn.min_nodes[i].top_to;
        vector<pair<int,int>> to_nbrs;
        if (edgeptr!=nullptr) {
            do {
                vector<pair<int,int>>::iterator it_find = find_if(to_nbrs.begin(),to_nbrs.end(), \
                    [edgeptr] (pair<int,int> elem) { return elem.first==edgeptr->from_node->min_id; } );
                if (it_find != to_nbrs.end()) {
                    int idx = distance(to_nbrs.begin(),it_find);
// for now let's just delete the newfound edge
//                    ktn.ts_edges[to_nbrs[idx].second].w += ktn.ts_edges[edgeptr->ts_pos].w;
//                    ktn.ts_edges[to_nbrs[idx].second+1].w += ktn.ts_edges[edgeptr->ts_pos+1].w; // ?
                    // what about the other edges?
                    ktn.del_spec_to_edge(edgeptr->to_node->min_id-1,edgeptr->ts_id);
                    ktn.del_spec_from_edge(edgeptr->to_node->min_id-1,edgeptr->ts_id);
                    ktn.del_spec_to_edge(edgeptr->from_node->min_id-1,edgeptr->ts_id);
                    ktn.del_spec_from_edge(edgeptr->from_node->min_id-1,edgeptr->ts_id);
                    ktn.n_edges = ktn.n_edges-2;
                } else { to_nbrs.emplace_back(make_pair(edgeptr->from_node->min_id,edgeptr->ts_pos));
                }
                edgeptr = edgeptr->next_to;
            } while (edgeptr!=nullptr);
        }
    }
    cout << "Finished checking for duplicate edges" << endl;
    // calculate in-degree and out-degree for nodes
    for (int i=0;i<nmin;i++) {
        double in_deg = Network::calc_deg_inout(ktn,i,1); ktn.min_nodes[i].deg_in = in_deg;
        double out_deg = Network::calc_deg_inout(ktn,i,2); ktn.min_nodes[i].deg_out = out_deg;
    }
}
