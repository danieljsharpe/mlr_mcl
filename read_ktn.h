// header file included in mlr_mcl.cpp

#include <vector>

using namespace std;

class Read_ktn {
    public:
    static vector<pair<int,int>> read_ts_conns(int);
    static vector<double> read_ts_weights(int);
};
