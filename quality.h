/* Functions for calculating metrics to assess the quality of a clustering and writing output */

#ifndef __QUALITY_H_INCLUDED__
#define __QUALITY_H_INCLUDED__

#include "ktn.h"

using namespace std;

class Quality_clust {

    public:

    static double calc_modularity(const Network&);
    static double calc_avgncut(const Network&);
    static void find_intercomm_edges(const Network&);
    static void write_comms(const Network&);
};

#endif
