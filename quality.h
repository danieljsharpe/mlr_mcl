/* Functions for calculating metrics to assess the quality of a clustering and writing output */

#ifndef __QUALITY_H_INCLUDED__
#define __QUALITY_H_INCLUDED__

#include <queue>
#include <utility>
#include <algorithm>
#include <numeric>
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
    static void check_if_disjoint(Network&,vector<int>,vector<int>);
    static void post_processing(Network&,int);
    static void post_processing_rates(Network&,int,int);
    static void set_new_comm_ids(Network&,vector<int>,int);
    static pair<vector<int>,vector<double>> ktn_bfs(Network&,int);
};

/* function to return a queued entry when popping from a queue */
template<class Q>
typename Q::value_type pop_from_queue(Q& q) {
    auto res = q.front();
    q.pop();
    return res;
};

/* function to get a vector of the ordered indices corresponding to sorted elems of vector */
template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {
  vector<size_t> idx_vec(v.size());
  iota(idx_vec.begin(), idx_vec.end(), 0);
  stable_sort(idx_vec.begin(), idx_vec.end(),
       [&v](size_t i1, size_t i2) { return v[i1]>v[i2]; }); // descending order
  return idx_vec;
}

#endif
