// header file included in mlr_mcl.cpp

#ifndef __READ_KTN_H_INCLUDED__
#define __READ_KTN_H_INCLUDED__

#include <vector>

using namespace std;

class Read_ktn {
    public:
    static vector<pair<int,int>> read_ts_conns(int);
    static vector<double> read_ts_weights(int);
};

#endif
