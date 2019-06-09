#include "quality.h"
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

/* Calculate the modularity objective function, Q. Append value to file "mod.dat" */
double Quality_clust::calc_modularity(const Network &ktn) {
    long double Q = 0.;
    long double m = 0.; // tot. weight of edges in the network
    for (int i=0;i<ktn.tot_edges;i++) {
        if (ktn.ts_edges[i].deadts) continue;
        m = log(exp(m) + exp(ktn.ts_edges[i].w));
    }
    m = exp(m);
    for (int i=0;i<ktn.tot_edges;i++) {
        if (ktn.ts_edges[i].deadts) continue;
        if (ktn.ts_edges[i].to_node->comm_id == ktn.ts_edges[i].from_node->comm_id) {
            Q += exp(ktn.ts_edges[i].w) - ((1./m)*\
                 (exp(ktn.ts_edges[i].to_node->deg_in)*exp(ktn.ts_edges[i].from_node->deg_out))); }
    }
    Q *= 1./m;
    ofstream mod_f;
    mod_f.open("mod.dat",ofstream::app);
    mod_f << Q << endl;
    mod_f.close();
    return Q;
}

/* Calculate the average normalised cut objective function, Ncut. Requires an additional file
   "stat_prob.dat" containing the stationary probabilities corresponding to the transition network.
   Append value to file "ncut.dat" */
double Quality_clust::calc_avgncut(const Network &ktn) {
    double Ncut = 0.;
    return Ncut;
}

/* Write file "sce_edge.dat" of 0/1 values representing whether an edge connects two communities.
   Average over many runs to determine probabilities and hence the clustering entropy */
void Quality_clust::find_intercomm_edges(Network &ktn) {

    ofstream sce_edge_f, sce_node_f;
    sce_edge_f.open("sce_edge.dat",ofstream::trunc);
    for (int i=0;i<ktn.tot_edges;i++) {
        if (i%2==1) continue; // edges are bidirected, skip second entry
        if (ktn.ts_edges[i].deadts) { // dead TS is indicated by -1 flag in output files
            sce_edge_f << -1 << endl; continue; }
        if (ktn.ts_edges[i].to_node->comm_id != ktn.ts_edges[i].from_node->comm_id) {
            sce_edge_f << 1 << endl;
            ktn.ts_edges[i].to_node->atboundary = true; ktn.ts_edges[i].from_node->atboundary = true;
        } else { sce_edge_f << 0 << endl; }
    }
    sce_edge_f.close();
    sce_node_f.open("sce_node.dat",ofstream::trunc);
    for (int i=0;i<ktn.tot_nodes;i++) {
        if (ktn.min_nodes[i].atboundary) { sce_node_f << 1 << endl;
        } else { sce_node_f << 0 << endl; }
    }
    sce_node_f.close();
}

/* Write communities to which nodes belong to a file "communities.dat" and write the
   attractor nodes to a file "attractors.dat" */
void Quality_clust::write_comms(const Network &ktn) {

    ofstream comms_f, attractors_f;
    comms_f.open("communities.dat",ofstream::trunc); attractors_f.open("attractors.dat",ofstream::trunc);
    for (int i=0;i<ktn.tot_nodes;i++) {
        if (ktn.min_nodes[i].comm_id==-1) {
            cout << "Error: node " << i << " has not been assigned to a community" << endl;
            throw Network::Ktn_exception(); }
        comms_f << ktn.min_nodes[i].comm_id << endl;
        if (ktn.min_nodes[i].attractor) attractors_f << ktn.min_nodes[i].min_id << endl;
    }
    comms_f.close(); attractors_f.close();
}
