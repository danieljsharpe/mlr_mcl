#include "quality.h"
#include <iostream>
#include <fstream>

using namespace std;

/* Calculate the modularity objective function */
double Quality_clust::calc_modularity(const Network &ktn) {

    return 1.;
}

/* Calculate the average normalised cut objective function. Requires an additional file
   "stat_prob.dat" containing the stationary probabilities corresponding to the transition network */
double Quality_clust::calc_avgncut(const Network &ktn) {

    return 1.;
}

/* Write a file of 0/1 values representing whether an edge connects two communities.
   Average over many runs to determine probabilities and hence the clustering entropy */
void Quality_clust::find_intercomm_edges(const Network &ktn) {

}

/* Write communities to which nodes belong to a file "communities.dat" and write the
   attractor nodes to a file "attractors.dat" */
void Quality_clust::write_comms(const Network &ktn) {

    ofstream comms_f, attractors_f;
    comms_f.open("communities.dat"); attractors_f.open("attractors.dat");
    for (int i=0;i<ktn.tot_nodes;i++) {
        if (ktn.min_nodes[i].comm_id==-1) {
            cout << "Error: node " << i << " has not been assigned to a community" << endl;
            throw Network::Ktn_exception(); }
        comms_f << ktn.min_nodes[i].comm_id << endl;
        if (ktn.min_nodes[i].attractor) attractors_f << ktn.min_nodes[i].min_id << endl;
    }
    comms_f.close(); attractors_f.close();
}
