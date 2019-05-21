#include "ktn.h"
#include <iostream>

using namespace std;

/* loop over TO or FROM edges downstream from edgeptr in linked list (used in debug tests) */
void print_edgeptr_info(Edge *edgeptr, int opt) {
    do {
        cout << "ts_id: " << edgeptr->ts_id << " w " << edgeptr->w << " FROM " << edgeptr->from_node->min_id << \
                " TO " << edgeptr->to_node->min_id << endl;
        if (opt==1) { edgeptr = edgeptr->next_to; } // scan TO nodes
        else if (opt==2) { edgeptr = edgeptr->next_from; } // scan FROM nodes
    } while (edgeptr != nullptr);
}

/* test merging of nodes */
void test_merge(Network &ktn, int n1, int n2) {
    Edge *edgeptr;
    cout << "\n\n*******\n\nmerging nodes " << n1 << " and " << n2 << "..." << endl;
    cout << "\noriginally edges pointing TO node " << n1 << ":" << endl;
    edgeptr = ktn.min_nodes[n1-1].top_to;
    print_edgeptr_info(edgeptr,1);
    cout << "\noriginally edges pointing FROM node " << n1 << ":" << endl;
    edgeptr = ktn.min_nodes[n1-1].top_from;
    print_edgeptr_info(edgeptr,2);
    cout << "\noriginally edges pointing TO node " << n2 << ":" << endl;
    edgeptr = ktn.min_nodes[n2-1].top_to;
    print_edgeptr_info(edgeptr,1);
    cout << "\noriginally edges pointing FROM node " << n2 << ":" << endl;
    edgeptr = ktn.min_nodes[n2-1].top_from;
    print_edgeptr_info(edgeptr,2);
    cout << "\nnow merge nodes " << n1 << " & " << n2 << "..." << endl;
    ktn.merge_nodes(n1-1,n2-1);
    cout << "\nnew edges pointing TO node " << n1 << ":" << endl;
    edgeptr = ktn.min_nodes[n1-1].top_to;
    print_edgeptr_info(edgeptr,1);
    cout << "\nnode " << n2 << " has no TO edges? " << (ktn.min_nodes[n2-1].top_to==nullptr) << \
            " node " << n2 << " is deleted? " << ktn.min_nodes[n2-1].deleted << endl;
    cout << "\n new edges pointing FROM node " << n1 << ":" << endl;
    edgeptr = ktn.min_nodes[n1-1].top_from;
    print_edgeptr_info(edgeptr,2);
    cout << "\nnode " << n2 << " has no FROM edges? " << (ktn.min_nodes[n2-1].top_from==nullptr) << endl;
}

/* debug tests for implementation of data structure */
void run_debug_tests(Network &ktn) {

    // TESTS
    // trivial tests
    cout << ktn.min_nodes[3].min_id << " (should be 4)" << endl;
    cout << ktn.ts_edges[2*500].from_node->min_id << " " << ktn.ts_edges[2*500].to_node->min_id << endl;
    cout << ktn.ts_edges[(2*500)+1].from_node->min_id << " " << ktn.ts_edges[(2*500)+1].to_node->min_id << endl;
    // loop over TO nodes
    cout << "iterate over all edges pointing TO node 1..." << endl;
    int i=0;
    Node *nodeptr; Edge *edgeptr;
    edgeptr = ktn.min_nodes[i].top_to;
    print_edgeptr_info(edgeptr,1);
    // loop over FROM nodes
    cout << "\niterate over all edges pointing FROM node " << i+1 << "..." << endl;
    edgeptr = ktn.min_nodes[i].top_from;
    print_edgeptr_info(edgeptr,2);
    // delete a single TO edge
    int j=252;
    cout << "\ndeleting TO edge with ts_id " << j << " from stack for node " << i+1 << "..." << endl;
    ktn.del_spec_to_edge(i,j);
    edgeptr = ktn.min_nodes[i].top_to;
    print_edgeptr_info(edgeptr,1);
    // update an edge so that it now is TO a new node
    j=3000;
    edgeptr = &ktn.ts_edges[j];
    cout << "\nupdating edge with ts_pos " << edgeptr->ts_pos << " so that it points TO node " << i+1 << endl;
    cout << "edge originally points FROM " << edgeptr->from_node->min_id << " TO " << \
            edgeptr->to_node->min_id << " ts_id: " << edgeptr->ts_id << endl;
    cout << "TO edges for the original TO node:" << endl;
    edgeptr = ktn.min_nodes[edgeptr->to_node->min_id-1].top_to;
    int old_to_min_id = edgeptr->to_node->min_id;
    print_edgeptr_info(edgeptr,1);
    cout << "now update TO edge:" << endl;
    ktn.update_to_edge(i,j);
    cout << "new TO edges for the original TO node:" << endl;
    edgeptr = ktn.min_nodes[old_to_min_id-1].top_to;
    print_edgeptr_info(edgeptr,1);
    cout << "new TO edges for the new TO node:" << endl;
    edgeptr = ktn.min_nodes[i].top_to;
    print_edgeptr_info(edgeptr,1);
    cout << "now specifically delete the top TO edge for the original TO node:" << endl;
    edgeptr = ktn.min_nodes[old_to_min_id-1].top_to;
    ktn.del_spec_to_edge(old_to_min_id-1,edgeptr->ts_id);
    cout << "  is top TO for this node nullptr? " << (ktn.min_nodes[old_to_min_id-1].top_to==nullptr) << endl;
    // delete all TO edges
    cout << "\nnow deleting all edges pointing TO node " << i+1 << "..." << endl;
    edgeptr = ktn.min_nodes[i].top_to;
    Edge **edgeptrptr = &ktn.min_nodes[i].top_to;
    cout << "edgeptrptr initially points to: " << (*edgeptrptr)->ts_id << endl;
    do {
        cout << "deleting edge, ts_id " << edgeptr->ts_id << endl;
        ktn.del_to_edge(i);
        edgeptr = edgeptr->next_to;
        if ((*edgeptrptr) != nullptr) {
            cout << "edgeptrptr now points to: " << (*edgeptrptr)->ts_id << endl;
        } else {
            cout << "edgeptrptr now points to nullptr" << endl;
        }
    } while (edgeptr != nullptr);
    // try to delete a TO edge now that such an edge no longer exists
    cout << "\ntry to delete an edge pointing TO node 1 (no such edge exists now)..." << endl;
    try {
        ktn.del_to_edge(i);
    } catch (Network::Ktn_exception& ktn_exc) {
        cout << "I caught a KTN_exception!" << endl;
    }

    // test merging of nodes
    test_merge(ktn,6,100); // merge nodes with min ID's 6 and 100
}
