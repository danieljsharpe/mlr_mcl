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

    cout << ">>>>> reading communities from file" << endl;
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
    vector<int> comm_sizes(ktn.n_comms);
    cout << ">>>>> no. of communities: " << n_comms << endl;
    for (auto &node: ktn.min_nodes) comm_sizes[node.comm_id]++;
    cout << ">>>>> reading attractors from file" << endl;
    vector<int> attractors = Read_ktn::read_single_col<int>(ktn.n_comms,"attractors.dat");
    for (auto att: attractors) {
        ktn.min_nodes[att-1].attractor=true; }
    check_if_disjoint(ktn,comm_sizes,attractors);
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

/* function checks that each community is a fully connected component */
void Quality_clust::check_if_disjoint(Network &ktn, vector<int> comm_sizes, vector<int> attractors) {
    // calculate community sizes and compare against results of BFS starting from attractor node
    for (int i=0;i<ktn.n_comms;i++) {
        pair<vector<int>,vector<double>> node_comm_masks = ktn_bfs(ktn,attractors[i]);
        vector<int> node_disc=node_comm_masks.first;
        int comm_size=0, j=0;
        for (int x: node_disc) {
            if (x && ktn.min_nodes[j].comm_id==i) comm_size++; j++; }
        cout << "  comm: " << i << "  comm size from file: " << comm_sizes[i] \
             << "  comm size from BFS: " << comm_size << endl;
        if (comm_sizes[i]!=comm_size) { // comm size from file does not equal comm size calculated from BFS
            cout << "    nodes not found in community from BFS: " << endl;
            int j=0;
            for (int x: node_disc) {
                if (!x && ktn.min_nodes[j].comm_id==i) cout << "  " << j+1; j++; }
            cout << endl; } // quack
//            throw Network::Ktn_exception(); }
    }
}

/* Communities below the size threshold are merged with the neighbouring community to which the community is connected
   by the fastest rate, unless the size of the neighbouring community exceeds the cap */
void Quality_clust::post_processing_rates(Network &ktn, int min_sz, int cap_sz) {

    cout << ">>>>> post-processing network to check rates to neighbours for communities of < " << min_sz << " nodes" << endl;
    int n_comms = ktn.n_comms; int merge_comm;
    vector<int> comm_sizes(n_comms);
    for (int i=0;i<ktn.tot_nodes;i++) {
        comm_sizes[ktn.min_nodes[i].comm_id]++; }
    for (int active_comm=0;active_comm<ktn.n_comms;active_comm++) {
        if (comm_sizes[active_comm]>=min_sz) continue; // community already exceeds minimum size
        int rep_node; int j=0;
        for (auto &node: ktn.min_nodes) {
            if (node.comm_id==active_comm) { rep_node=j+1; break; }
            j++;
        }
        pair<vector<int>,vector<double>> node_comm_masks = ktn_bfs(ktn,rep_node);
        vector<double> rates_to_comms = node_comm_masks.second;
        vector<long unsigned int> idx_comm_trans = sort_indexes<double>(rates_to_comms);
        merge_comm=-1; int idx;
        // ignore first elem in following loop, the fastest rate is a dummy value representing the current node
        for (j=1;j<ktn.n_comms;j++) {
            idx=idx_comm_trans[j];
            if (comm_sizes[idx]>cap_sz) { continue;
            } else if (rates_to_comms[idx]==-numeric_limits<double>::infinity()) { break;
            } else { merge_comm=idx; break; }
        }
        if (merge_comm==-1) continue; // there is no suitable neighbouring community for merging
        cout << "  merging comm " << active_comm << "  of size " << comm_sizes[active_comm] << "    with    comm " \
             << merge_comm << "  of size " << comm_sizes[merge_comm] << endl;
        comm_sizes[merge_comm] += comm_sizes[active_comm]; comm_sizes[active_comm] = 0;
        n_comms--;
        for (auto &node: ktn.min_nodes) { // set new communities and attractors in the ktn data structure
            if (node.comm_id!=active_comm) continue;
            node.comm_id = merge_comm;
            if (node.attractor) node.attractor = false; }
    }
    set_new_comm_ids(ktn,comm_sizes,n_comms);
}

/* The smallest communities are repeatedly merged with the smallest neighbouring communities, until the
   number of nodes in all communities exceeds min_sz */
void Quality_clust::post_processing(Network &ktn, int min_sz) {
    cout << ">>>>> post-processing network to subsume communities of < " << min_sz << " nodes" << endl;
    int n_comms = ktn.n_comms; // keep track of number of communities
    vector<int> comm_sizes(n_comms); // number of nodes in each community
    int active_comm; // ID of the community whose members are being searched for (current smallest community)
    int merge_comm; // ID of the community to which the active community is to be merged
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
        int rep_node=-1; // representative node of current smallest community
        int j=0;
        for (auto &node: ktn.min_nodes) {
            if (node.comm_id==active_comm) { rep_node=j+1; break; }
            j++;
        }
        pair<vector<int>,vector<double>> node_comm_masks = ktn_bfs(ktn,rep_node);
        vector<int> comm_disc(ktn.n_comms); // flags for neighbouring communities     node_comm_masks.second;
        for (int i=0;i<ktn.n_comms;i++) {
            if (node_comm_masks.second[i]>-numeric_limits<double>::infinity()) comm_disc[i]=1;
        }
        j=0; comm_disc[active_comm]=0;
        merge_comm=distance(begin(comm_disc),find_if(begin(comm_disc),end(comm_disc),[](int x) { return x != 0; }));
        for (auto comm_nbr: comm_disc) {
            if (comm_nbr && comm_sizes[j]<comm_sizes[merge_comm]) merge_comm=j; j++; }

        cout << "  merging comm " << active_comm << "  of size " << comm_sizes[active_comm] << "    with    comm " \
             << merge_comm << "  of size " << comm_sizes[merge_comm] << endl;
        comm_sizes[merge_comm] += comm_sizes[active_comm]; comm_sizes[active_comm] = 0;
        n_comms--;
        for (auto &node: ktn.min_nodes) { // set new communities and attractors in the ktn data structure
            if (node.comm_id!=active_comm) continue;
            node.comm_id = merge_comm;
            if (node.attractor) node.attractor = false; }
        active_comm = min_element(comm_sizes.begin(),comm_sizes.end(),get_smallest_nonzero) - comm_sizes.begin();
    }
    set_new_comm_ids(ktn,comm_sizes,n_comms);
}

/* trace community sizes and update the comm_id's of Node structures to be consistent with the new numbering */
void Quality_clust::set_new_comm_ids(Network &ktn, vector<int> comm_sizes, int n_comms) {
    cout << ">>>>> setting new community IDs" << endl;
    int n_new_comms=0;
    Node *nodeptr;
    unordered_map<int,int> comm_ids_map;
    comm_ids_map.reserve(ktn.n_comms);
    for (int i=0;i<ktn.n_comms;i++) {
        if (comm_sizes[i]!=0) { comm_ids_map.insert({i,n_new_comms}); n_new_comms++; }
    }
    for (int i=0;i<ktn.n_nodes;i++) {
        nodeptr = &ktn.min_nodes[i];
        int new_comm_id = comm_ids_map.at(nodeptr->comm_id);
        nodeptr->comm_id = new_comm_id;
    }
    ktn.n_comms=n_comms;
    cout << ">>>>> no. of communities after merging: " << ktn.n_comms << endl;
}

/* use breadth first search, initiated from a particular node (eg the attractor of the community),
   to find all members of a given community */
pair<vector<int>,vector<double>> Quality_clust::ktn_bfs(Network &ktn, int att_min_id) {

    int active_comm = ktn.min_nodes[att_min_id-1].comm_id;
    queue<int> nbr_queue; // queue of nodes of the active community to visit
    vector<int> node_disc(ktn.tot_nodes); // flags to indicate BFS has discovered a node
    // log rates of transitions to neighbouring communities. Used as flags to indicate BFS has discovered a neighbouring community
    vector<double> comm_disc(ktn.n_comms,-numeric_limits<double>::infinity());
    node_disc[att_min_id-1]=1;
    comm_disc[active_comm]=numeric_limits<double>::infinity(); // dummy value indicates that this is the current community
    nbr_queue.push(att_min_id); // initial node to start BFS procedure
    while (!nbr_queue.empty()) {
        int node_id = pop_from_queue(nbr_queue);
        Node *nodeptr = &ktn.min_nodes[node_id-1];
        Edge *edgeptr = nodeptr->top_from;
        while (edgeptr != nullptr) {
            if (edgeptr->to_node->comm_id!=active_comm) {
                comm_disc[edgeptr->to_node->comm_id] = log(exp(comm_disc[edgeptr->to_node->comm_id])+exp(edgeptr->w));
                node_disc[edgeptr->to_node->min_id-1] = 1;
            } else if (edgeptr->to_node->comm_id == active_comm && !node_disc[edgeptr->to_node->min_id-1]) {
                node_disc[edgeptr->to_node->min_id-1] = 1;
                nbr_queue.push(edgeptr->to_node->min_id);
            }
        edgeptr = edgeptr->next_from;
        }
    }
    return make_pair(node_disc,comm_disc);
}
