#include "quality.h"
#include "read_ktn.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <limits>

using namespace std;

/* Calculate the modularity objective function, Q. Append value to file "modularity.dat".
   This function processes the unclustered network (cannot use this function after merging edge weights -
   hence do this calculation in separate execution and after any desired post-processing) */
double Quality_clust::calc_modularity(const Network &ktn) {

    long double Q = 0.;
    for (int i=0;i<ktn.tot_nodes;i++) {
        for (int j=0;j<ktn.tot_nodes;j++) {
            double w_ij = -numeric_limits<double>::infinity();
            Edge *edgeptr = ktn.min_nodes[i].top_from;
            if (i==j) goto increment_mod; // transition rate matrix has no self-loops
            if (ktn.min_nodes[i].comm_id != ktn.min_nodes[j].comm_id) continue;
            if (edgeptr==nullptr) {
                cout << "Error: node " << i+1 << " has no neighbours" << endl;
                throw Network::Ktn_exception(); }
            do {
               if (edgeptr->to_node->min_id-1==j) {
                   w_ij = edgeptr->w_r;
                   break; }
               edgeptr = edgeptr->next_from;
            } while (edgeptr!=nullptr);
            increment_mod:
                Q += exp(w_ij) - ((1./exp(ktn.tot_in_deg))*\
                     (exp(ktn.min_nodes[i].deg_in + ktn.min_nodes[j].deg_out)));
            }
    }
    Q *= 1./exp(ktn.tot_in_deg);
    ofstream mod_f;
    mod_f.open("modularity.dat",ofstream::app);
    mod_f << Q << endl;
    mod_f.close();
    return Q;
}

/* Calculate the average normalised cut objective function, Ncut. Requires that stationary occupation
   probabilities of minima have been read in to the ktn from the file "stat_prob.dat".
   Append value to file "avgncut.dat" */
double Quality_clust::calc_avgncut(const Network &ktn) {

    double Ncut = 0.;
    vector<double> stat_prob_comms(ktn.n_comms,-numeric_limits<double>::infinity());
    vector<double> Ncut_c(ktn.n_comms,0.); // Ncut values for each cluster
    for (int i=0;i<ktn.tot_nodes;i++) {
        stat_prob_comms[ktn.min_nodes[i].comm_id] = \
            log(exp(stat_prob_comms[ktn.min_nodes[i].comm_id]) + exp(ktn.min_nodes[i].peq));
    }
    for (int i=0;i<ktn.tot_edges;i++) {
        if (ktn.ts_edges[i].to_node->comm_id != ktn.ts_edges[i].from_node->comm_id) {
            Ncut_c[ktn.ts_edges[i].from_node->comm_id] += \
                (exp(ktn.ts_edges[i].from_node->peq - ktn.ts_edges[i].from_node->deg_out + \
                     ktn.ts_edges[i].w_r)/exp(stat_prob_comms[ktn.ts_edges[i].from_node->comm_id]));
            Ncut_c[ktn.ts_edges[i].to_node->comm_id] += \
                (exp(ktn.ts_edges[i].from_node->peq - ktn.ts_edges[i].from_node->deg_out + \
                     ktn.ts_edges[i].w_r)/(1.-exp(stat_prob_comms[ktn.ts_edges[i].to_node->comm_id])));
        }
    }
    for (int i=0;i<ktn.n_comms;i++) {
        Ncut += Ncut_c[i];
        cout << "comm i: " << i << "    avg n_cut: " << Ncut_c[i] << endl;
    }
    Ncut *= (1./double(ktn.n_comms));
    return Ncut;
}

/* Calculate the conductance objective function. Append value to file "conductance.dat" */
double Quality_clust::calc_conductance(const Network &ktn) {

    double conduc = -numeric_limits<double>::infinity();
    return conduc;
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

/* Read communities from output file (used when executing for post-processing) */
void Quality_clust::read_comms(Network &ktn) {

    vector<int> communities = Read_ktn::read_single_col<int>(ktn.tot_nodes,"communities.dat");
    int i=0, n_comms=0;
    for (const auto& comm: communities) {
        if (comm<0) {
            cout << "Error: node " << i+1 << " has not been assigned a valid community" << endl;
            throw Network::Ktn_exception(); }
        ktn.min_nodes[i].comm_id = comm; i++;
        if (comm+1>n_comms) n_comms = comm+1;
    }
    ktn.n_comms = n_comms;
}

/* Write communities to which nodes belong to a file "communities.dat" and write the
   attractor nodes to a file "attractors.dat" */
void Quality_clust::write_comms(const Network &ktn) {

    ofstream comms_f, attractors_f;
    comms_f.open("communities.dat",ofstream::trunc); attractors_f.open("attractors.dat",ofstream::trunc);
    for (int i=0;i<ktn.tot_nodes;i++) {
        if (ktn.min_nodes[i].comm_id<0) {
            cout << "Error: node " << i+1 << " has not been assigned to a community" << endl;
            throw Network::Ktn_exception(); }
        comms_f << ktn.min_nodes[i].comm_id << endl;
        if (ktn.min_nodes[i].attractor) attractors_f << ktn.min_nodes[i].min_id << endl;
    }
    comms_f.close(); attractors_f.close();
}

void Quality_clust::post_processing(Network &ktn, int min_sz) {
    cout << ">>>>> post-processing network to subsume communities of <" << min_sz << " nodes" << endl;
}
