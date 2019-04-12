/*
Custom data structure for representing and reading in a ktn
*/

#include <exception>
#include <vector>

using namespace std;

struct Node;

struct Edge {
    int ts_id;
    int ts_pos; // position of the TS in the ts_edges array
    double w;
    bool deadts = false;
    Node *to_node;
    Node *from_node;
    Edge *next_to;
    Edge *next_from;
};

struct Node {
    int min_id;
    bool hem_flag; // flag for use in heavy edge matching routine
    bool deleted = false; // indicates node has been "deleted" from the network
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
    void merge_nodes(int,int);
    static void setup_network(Network&,int,int,vector<pair<int,int>>,vector<double>);
    vector<Node> min_nodes;
    vector<Edge> ts_edges;

    struct Ktn_exception {
        const char * what () const throw () { return "KTN Exception"; }
    };

    int n_nodes; int n_edges;
};
