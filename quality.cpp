#include "quality.h"
#include "read_ktn.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <limits>
#include <algorithm>
#include <functional>
#include <unordered_map>

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
        if (ktn.ts_edges[i].deadts) continue;
        if (ktn.ts_edges[i].to_node->comm_id != ktn.ts_edges[i].from_node->comm_id) {
            Ncut_c[ktn.ts_edges[i].from_node->comm_id] += \
                (exp(ktn.ts_edges[i].from_node->peq - ktn.ts_edges[i].from_node->deg_out + \
                     ktn.ts_edges[i].w_r)/exp(stat_prob_comms[ktn.ts_edges[i].from_node->comm_id]));
        }
    }
    for (int i=0;i<ktn.n_comms;i++) {
        Ncut += Ncut_c[i];
    }
    Ncut *= (1./double(ktn.n_comms));
    ofstream avgncut_f;
    avgncut_f.open("avgncut.dat",ofstream::app);
    avgncut_f << Ncut << endl;
    avgncut_f.close();
    return Ncut;
}

/* Calculate the conductance objective function. Append value to file "conductance.dat" */
double Quality_clust::calc_conductance(const Network &ktn) {

    vector<double> deg_out_comms(ktn.n_comms,-numeric_limits<double>::infinity());
    vector<double> conduc_c(ktn.n_comms,-numeric_limits<double>::infinity()); // conductance values for each community
    for (int i=0;i<ktn.tot_nodes;i++) {
        deg_out_comms[ktn.min_nodes[i].comm_id] = \
            log(exp(deg_out_comms[ktn.min_nodes[i].comm_id]) + exp(ktn.min_nodes[i].deg_out));
    }
    for (int i=0;i<ktn.tot_edges;i++) {
        if (ktn.ts_edges[i].deadts) continue;
        if (ktn.ts_edges[i].to_node->comm_id != ktn.ts_edges[i].from_node->comm_id) {
            conduc_c[ktn.ts_edges[i].from_node->comm_id] = \
                log(exp(conduc_c[ktn.ts_edges[i].from_node->comm_id]) + exp(ktn.ts_edges[i].w_r));
        }
    }
    // lambda function to find the minimum volume - either of the community or its complement
    auto get_conduc_denom = [](double deg_out_comm, double deg_out_tot) {
        // (log) out-degree of set that is complement of the community
        double deg_out_comp = log(exp(deg_out_tot)-exp(deg_out_comm));
        if (deg_out_comm<deg_out_comp) { return exp(deg_out_comm); }
        else { return exp(deg_out_comp); }
    };
    int c_i=0; double conduc=0.;
    using funcptr = double(*)(double,double);
    funcptr fptr = get_conduc_denom;
    for (auto &phi_c: conduc_c) {
        phi_c = exp(phi_c) / get_conduc_denom(deg_out_comms[c_i],ktn.tot_in_deg);
        conduc += phi_c;
        c_i++;
    }
    conduc *= (1./double(ktn.n_comms));
    ofstream conduc_f;
    conduc_f.open("conductance.dat",ofstream::app);
    conduc_f << conduc << endl;
    conduc_f.close();
    return conduc;
}

/* Write file "sce_edge.dat" of 0/1 values representing whether an edge connects two communities.
   Average over many runs to determine probabilities and hence the clustering entropy.
   Write file "sce_node.dat" representing whether a node is at the boundary between two communities. */
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

/* Read communities and attractors from previous output files (used when executing for post-processing) */
void Quality_clust::read_comms(Network &ktn) {

    vector<int> communities = Read_ktn::read_single_col<int>(ktn.tot_nodes,"communities.dat");
    int i=0, n_comms=0;
    for (auto comm: communities) {
        if (comm<0) {
            cout << "Error: node " << i+1 << " has not been assigned a valid community" << endl;
            throw Network::Ktn_exception(); }
        ktn.min_nodes[i].comm_id = comm; i++;
        if (comm+1>n_comms) n_comms = comm+1;
    }
    ktn.n_comms = n_comms;
    vector<int> attractors = Read_ktn::read_single_col<int>(ktn.n_comms,"attractors.dat");
    for (auto att: attractors) {
        ktn.min_nodes[att-1].attractor=true; }
}

/* Write communities to which nodes belong to a file "communities.dat" and write the
   attractor nodes to a file "attractors.dat" */
void Quality_clust::write_comms(const Network &ktn, int opt) {

    ofstream comms_f, attractors_f;
    if (opt==0) {
        comms_f.open("communities.dat",ofstream::trunc); attractors_f.open("attractors.dat",ofstream::trunc);
    } else if (opt==1) {
        comms_f.open("communities_new.dat",ofstream::trunc);
        attractors_f.open("attractors_new.dat",ofstream::trunc);
    }
    for (int i=0;i<ktn.tot_nodes;i++) {
        if (ktn.min_nodes[i].comm_id<0) {
            cout << "Error: node " << i+1 << " has not been assigned to a community" << endl;
            throw Network::Ktn_exception(); }
        comms_f << ktn.min_nodes[i].comm_id << endl;
        if (ktn.min_nodes[i].attractor) attractors_f << ktn.min_nodes[i].min_id << endl;
    }
    comms_f.close(); attractors_f.close();
}

/* The smallest communities are repeatedly merged with the smallest neighbouring communities, until the
   number of nodes in all communities exceeds min_sz */
void Quality_clust::post_processing(Network &ktn, int min_sz) {
    cout << ">>>>> post-processing network to subsume communities of < " << min_sz << " nodes" << endl;
    int n_comms = ktn.n_comms;
    vector<double> comm_sizes(n_comms); // number of nodes in each community
    int active_comm; // ID of the community whose members are being searched for
    int merge_comm; // ID of the community to which the active community is to be merged
    int n_disc;
    Node *nodeptr; Edge *edgeptr;
    for (int i=0;i<ktn.tot_nodes;i++) {
        comm_sizes[ktn.min_nodes[i].comm_id]++; }
    // lambda comparator returning bool to indicate if a<b, if a is nonzero
    function<bool(int a, int b)> get_smallest_nonzero = [](int a, int b) {
        if ((a==0) && (b!=0)) { return false;
        } else if ((a!=0) && (b==0)) { return true;
        } else if (a<b) { return true;
        } else { return false; }
    };
    active_comm = min_element(comm_sizes.begin(),comm_sizes.end(),get_smallest_nonzero) - comm_sizes.begin();
    if (comm_sizes[active_comm] >= min_sz) {
        cout << "    all communities already exceed the size threshold" << endl; return; }
    while (comm_sizes[active_comm] < min_sz) {
        // breadth-first search (BFS) to find all nodes of the active community
        queue<int> nbr_queue; // queue of nodes of the active community to visit
        vector<int> node_disc(ktn.tot_nodes); // flags to indicate BFS has discovered a node
        vector<int> comm_disc(n_comms); // flags to indicate BFS has discovered a (neighbouring) community
        vector<int> nbr_comm, nbr_nodes; // list of neighbouring communities and nodes, respectively
        n_disc=0;
        for (int i=0;i<ktn.tot_nodes;i++) {
            if (ktn.min_nodes[i].comm_id==active_comm) {
                node_disc[ktn.min_nodes[i].min_id-1] = 1; n_disc++;
                nbr_nodes.emplace_back(ktn.min_nodes[i].min_id);
                nbr_queue.push(ktn.min_nodes[i].min_id); break; }
        }
        if (n_disc==0) {
            cout << "Error: selected community " << active_comm << " has no members" << endl;
            throw Network::Ktn_exception(); }
        while (!nbr_queue.empty()) {
            int node_id = pop_from_queue(nbr_queue);
            nodeptr = &ktn.min_nodes[node_id-1];
            edgeptr = nodeptr->top_from;
            while (edgeptr != nullptr) {
                if ((edgeptr->to_node->comm_id != active_comm) && (!comm_disc[edgeptr->to_node->comm_id])) {
                    comm_disc[edgeptr->to_node->comm_id] = 1;
                    nbr_comm.emplace_back(edgeptr->to_node->comm_id);
                } else if ((edgeptr->to_node->comm_id == active_comm) && (!node_disc[edgeptr->to_node->min_id-1])) {
                    node_disc[edgeptr->to_node->min_id-1] = 1; n_disc++;
                    nbr_nodes.emplace_back(edgeptr->to_node->min_id);
                    nbr_queue.push(edgeptr->to_node->min_id);
                }
            edgeptr = edgeptr->next_from;
            }
        }
        if (n_disc != comm_sizes[active_comm]) {
            cout << "Error: community " << active_comm << " is disjoint" << endl;
            throw Network::Ktn_exception(); }
        merge_comm = nbr_comm[0];
        for (auto comm: nbr_comm) {
            if (comm_sizes[comm] < comm_sizes[merge_comm]) merge_comm = comm; }
        cout << "  merging comm " << active_comm << "  size " << comm_sizes[active_comm] << "    with    comm " \
             << merge_comm << "  size " << comm_sizes[merge_comm] << endl;
        comm_sizes[merge_comm] += comm_sizes[active_comm]; comm_sizes[active_comm] = 0;
        ktn.n_comms--;
        for (auto node: nbr_nodes) { // set new communities and attractors in the ktn data structure
            ktn.min_nodes[node-1].comm_id = merge_comm;
            if (ktn.min_nodes[node-1].attractor) ktn.min_nodes[node-1].attractor = false; }
        active_comm = min_element(comm_sizes.begin(),comm_sizes.end(),get_smallest_nonzero) - comm_sizes.begin();
    }
    cout << ">>>>> no. of communities after merging: " << ktn.n_comms << endl;
    // trace comm_sizes and update comm_id's to be consistent
    int n_new_comms=0;
    unordered_map<int,int> comm_ids_map;
    comm_ids_map.reserve(ktn.n_comms);
    for (int i=0;i<n_comms;i++) {
        if (comm_sizes[i]!=0) { comm_ids_map.insert({i,n_new_comms}); n_new_comms++; }
    }
    for (int i=0;i<ktn.n_nodes;i++) {
        nodeptr = &ktn.min_nodes[i];
        int new_comm_id = comm_ids_map.at(nodeptr->comm_id);
        nodeptr->comm_id = new_comm_id;
    }
}
