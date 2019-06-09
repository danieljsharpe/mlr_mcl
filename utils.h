#ifndef __UTILS_H_INCLUDED__
#define __UTILS_H_INCLUDED__

#include "ktn.h"
#include <random>
#include <typeinfo>
#include <vector>

using namespace std;

/* Head file containing general utility functions NB template funcs have to be implemented
   inside the header so that the compiler can generate code for all specialisations of the code */

/* uniform random number generator given a range */
template <class T>
static T rand_unif(T xmin, T xmax, int seed) {
    static default_random_engine generator (seed);
    if (typeid(T) == typeid(double)) {
        uniform_real_distribution<double> distribution1(xmin,xmax);
        return distribution1(generator); }
    else if (typeid(T) == typeid(int)) {
        uniform_int_distribution<int> distribution2(xmin,xmax);
        return distribution2(generator); }
}

/* flatten a vector of vectors efficiently */
template <typename T>
vector<T> flatten(const vector<vector<T>>& vec) {
    size_t tot_sz = 0;
    for (const auto& sub: vec) // loop over subvectors
        tot_sz += sub.size();
    vector<T> vec_ret;
    vec_ret.reserve(tot_sz);
    for (const auto& sub: vec)
        vec_ret.insert(vec_ret.end(),sub.begin(),sub.end());
    return vec_ret;
}

void print_edgeptr_info(Edge *,int);
void run_debug_tests(Network&);

typedef pair<vector<pair<double,int>>,vector<int>> Csr_mtx; // matrix in CSR (or CSC) sparse format
void print_sparse_matrix(const Csr_mtx&);

#endif
