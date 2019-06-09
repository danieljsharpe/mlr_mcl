/*
Read in kinetic transition network information from files "ts_conns.dat", "ts_weights.dat" and "stat_prob.dat"
*/

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include "read_ktn.h"

using namespace std;

vector<pair<int,int>> Read_ktn::read_double_col(int nlines, string inp_fname) {

    string line;
    ifstream inp_f;
    inp_f.open(inp_fname);
    vector<pair<int,int>> entries_pairs(nlines);
    for (int i=0;i<nlines;i++) {
        getline(inp_f,line);
        istringstream min_pair(line);
        auto itt = istream_iterator<string>(min_pair);
        int min1 = stoi(*itt);
        itt++;
        int min2 = stoi(*itt);
        entries_pairs[i] = make_pair(min1,min2);
    }
    inp_f.close();
    return entries_pairs;
}

vector<double> Read_ktn::read_single_col(int nlines, string inp_fname) {

    string val_str;
    vector<double> entries(nlines);
    ifstream inp_f;
    inp_f.open(inp_fname);
    for (int i=0;i<nlines;i++) {
        getline(inp_f,val_str);
        entries[i] = stod(val_str);
    }
    inp_f.close();
    return entries;
}
