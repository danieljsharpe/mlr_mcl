// header file included in mlr_mcl.cpp

#ifndef __READ_KTN_H_INCLUDED__
#define __READ_KTN_H_INCLUDED__

#include <vector>
#include <string>

using namespace std;

class Read_ktn {
    public:
    static vector<pair<int,int>> read_double_col(int,string);
    static vector<double> read_single_col(int,string);
};

#endif
