/*
Read in kinetic transition network information from files "ts_conns.dat" and "ts_weights.dat"
*/

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include "read_ktn.h"

using namespace std;

vector<pair<int,int>> Read_ktn::read_ts_conns(int nts) {

    string line;
    ifstream ts_conns_f;
    ts_conns_f.open("ts_conns.dat");
    vector<pair<int,int>> ts_conns(nts);
    for (int i=0;i<nts;i++) {
        getline(ts_conns_f,line);
        istringstream min_pair(line);
        auto itt = istream_iterator<string>(min_pair);
        int min1 = stoi(*itt);
        itt++;
        int min2 = stoi(*itt);
        ts_conns[i] = make_pair(min1,min2);
    }
    ts_conns_f.close();
    return ts_conns;
}

vector<double> Read_ktn::read_ts_weights(int nts) {

    string ts_wt_str;
    vector<double> ts_weights(2*nts);
    ifstream ts_wts_f;
    ts_wts_f.open("ts_weights.dat");
    for (int i=0;i<2*nts;i++) {
        getline(ts_wts_f,ts_wt_str);
        ts_weights[i] = stod(ts_wt_str);
    }
    ts_wts_f.close();
    return ts_weights;
}
