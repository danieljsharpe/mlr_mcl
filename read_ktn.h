/*
Read in kinetic transition network information from files "ts_conns.dat", "ts_weights.dat" and "stat_prob.dat"
*/

#ifndef __READ_KTN_H_INCLUDED__
#define __READ_KTN_H_INCLUDED__

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <typeinfo>

using namespace std;

class Read_ktn {

    public:

    template <typename T> static vector<pair<T,T>> read_double_col(int,string);
    template <typename T> static vector<T> read_single_col(int,string);
};

template <typename T>
vector<pair<T,T>> Read_ktn::read_double_col(int nlines, string inp_fname) {

    string line;
    ifstream inp_f;
    inp_f.open(inp_fname);
    vector<pair<T,T>> entries_pairs(nlines);
    for (int i=0;i<nlines;i++) {
        T min1, min2;
        getline(inp_f,line);
        istringstream min_pair(line);
        auto itt = istream_iterator<string>(min_pair);
        if (typeid(T)==typeid(int)) { min1 = stoi(*itt);
        } else if (typeid(T)==typeid(double)) { min1 = stod(*itt); }
        itt++;
        if (typeid(T)==typeid(int)) { min2 = stoi(*itt);
        } else if (typeid(T)==typeid(double)) { min2 = stod(*itt); }
        entries_pairs[i] = make_pair(min1,min2);
    }
    inp_f.close();
    return entries_pairs;
}

template <typename T>
vector<T> Read_ktn::read_single_col(int nlines, string inp_fname) {

    string val_str;
    vector<T> entries(nlines);
    ifstream inp_f;
    inp_f.open(inp_fname);
    for (int i=0;i<nlines;i++) {
        getline(inp_f,val_str);
        if (typeid(T)==typeid(double)) { entries[i] = stod(val_str);
        } else if (typeid(T)==typeid(int)) { entries[i] = stoi(val_str); }
    }
    inp_f.close();
    return entries;
}

#endif
