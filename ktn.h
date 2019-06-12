/*
Custom data structure for representing and reading in a ktn
*/

#ifndef __KTN_H_INCLUDED__
#define __KTN_H_INCLUDED__

#include <exception>
#include <vector>

using namespace std;

struct Node;

struct Edge {
    int ts_id;
    int ts_pos; // position of the TS in the ts_edges array
    double w;
    double w_r; // original weight
    bool deadts = false; // indicates TS only linked to one minimum
    bool intercomm = false; // flag if edge connects nodes of two different communities
    Node *to_node;
    Node *from_node;
    Edge *next_to;
    Edge *next_from;
};

struct Node {
    int min_id;
    int comm_id = -1; // community ID (-1 indicates null value)
    bool hem_flag; // flag for use in heavy edge matching routine
    bool deleted = false; // indicates node has been "deleted" from the network
    bool attractor = false; // indicates node is an attractor (in MLR-MCL)
    bool atboundary = false; // indicates node has at least one inter-community edge
    double deg_in, deg_out; // (log) in-degree and out-degree
    double peq; // equilibrium (stationary) occupation probability
    Edge *top_to;
    Edge *top_from;
};

// structure containing the kinetic transition network
struct Network {

    public:

    Network(int,int);
    ~Network();
    void del_node(int);
    void add_to_edge(int,int);
    void add_from_edge(int,int);
    void del_to_edge(int);
    void del_from_edge(int);
    void del_spec_to_edge(int,int);
    void del_spec_from_edge(int,int);
    void update_to_edge(int,int);
    void update_from_edge(int,int);
    void merge_nodes(int,int);
    static double calc_deg_inout(const Network&,int,int);
    static void setup_network(Network&,int,int,const vector<pair<int,int>>,const vector<double>, \
        const vector<double>);
    vector<Node> min_nodes;
    vector<Edge> ts_edges;

    struct Ktn_exception {
        const char * what () const throw () { return "KTN Exception"; }
    };

    int n_nodes; int n_edges; int tot_nodes; int tot_edges; int n_dead=0;
    double tot_in_deg; // (log) total in-degree of nodes (N.B. = (log) total out-degree)
    int n_comms; // number of communities
};

#endif
