/* Functions for calculating metrics to assess the quality of a clustering and writing output */

#ifndef __QUALITY_H_INCLUDED__
#define __QUALITY_H_INCLUDED__

#include <queue>
#include "ktn.h"

using namespace std;

class Quality_clust {

    public:

    static double calc_modularity(const Network&);
    static double calc_avgncut(const Network&);
    static double calc_conductance(const Network&);
    static void find_intercomm_edges(Network&);
    static void read_comms(Network&);
    static void write_comms(const Network&,int);
    static void post_processing(Network&,int);
};

/* function to return a queued entry when popping from a queue */
template<class Q>
typename Q::value_type pop_from_queue(Q& q) {
    auto res = q.front();
    q.pop();
    return res;
};

#endif
