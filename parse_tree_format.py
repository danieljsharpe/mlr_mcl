'''
Python script to parse Infomap output in .tree format, ie:
<first level hierarchy idx>:<second level hierarchy index>:<third level hierarchy index>:... \
<steady state population of random walkers for node> "<new node name>" <old node name>

The 3-level hierarchy is:
major clusters / minor clusters / indiv nodes within minor cluster

The script then uses networkx package to plot the Infomap result

Execute with e.g.:
python parse_tree_format.py dolphin_network 40 1 2 1 2 0
This plots the 40 nodes of the network printed in dolphin_network.tree. Does not read a .flow file.
Finds 1st + 2nd level clusters. Plots individual nodes and colours them by second-level cluster IDs.

e.g.2.
python parse_tree_format.py dolphin_network 40 2 1 1 2 2
This time, read a .ftree format, read a flow file, find first + 2nd level clusters,
and use this info to plot & colour the second-level clusters.

Daniel J. Sharpe
Dec 2018
'''

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from collections import OrderedDict

class Parse_treefile(object):

    def __init__(self,rootfilename,n_nodes,treeformat,read_clust_level,readflow=True,plot_level=1):

        if treeformat==1: self.treefilename = rootfilename+".tree"
        elif treeformat==2: self.treefilename = rootfilename+".ftree" # also includes link weight info
        self.readflow = readflow
        self.rootfilename = rootfilename
        self.n_nodes = n_nodes
        self.read_clust_level = read_clust_level
        self.plot_level = plot_level
        self.parsefile()

    ''' Function to parse the .tree file '''
    def parsefile(self):

        self.weights = [] # weight of each indiv node
        self.level1clusts = OrderedDict([]) # list of nodes belonging to "level 1" clusters
        self.level2clusts = OrderedDict([]) # list of nodes belonging to "level 2" clusters
        self.level1clustflow = [] # tot. weight (internal flow) of each level 1 cluster
        self.level2clustflow = [] # tot. weight of each level 2 cluster
        self.node_clust1_id = [0]*self.n_nodes # list of "level 1" clusters to which each node belongs
        self.node_clust2_id = [0]*self.n_nodes # list of "level 2" clusters to which each node belongs
        current_clust1, current_clust2 = 0, 0
        if self.read_clust_level != 1:
            current_clusttop = 0 # we must keep track of top level cluster
            self.tot_mainclusters = 0
        self.tot_subclusters = 0
        cum_flow = 0. # cumulative total flow

        # parse .tree/.ftree file
        with open(self.treefilename,"r") as treef:
            for line in treef.readlines():
                line = line.split()
                if line[0] == "#": continue # comment line
                elif line[0] == "*Links": break # reached end of node info in ftree file
                hierarchy = line[0].split(":")
                # note: the no. of levels of hierarchy (>=2) can be variable even within a single .tree file
                try:
                    level1clust, level2clust = int(hierarchy[self.read_clust_level-1]), \
                                               int(hierarchy[self.read_clust_level])
                    if self.read_clust_level != 1: leveltopclust = int(hierarchy[0])
                except IndexError:
                    print "Error: trying to read a cluster level that doesn't exist for cluster ", hierarchy[0], \
                          "Try decreasing read_clust_level"
                    raise
                # NB need to account for the fact that we could have e.g. 1:1:1 and on next line 2:1:1, hence OR
                if level1clust != current_clust1 or \
                   (self.read_clust_level != 1 and leveltopclust != current_clusttop):
                        if self.read_clust_level == 1:
                            current_clust1 += 1
                            self.level1clusts[current_clust1] = []
                        else:
                            if leveltopclust != current_clusttop: current_clusttop += 1
                            current_clust1 = level1clust
                            self.tot_mainclusters += 1
                            self.level1clusts[self.tot_mainclusters] = []
                        self.level1clustflow.append(0.)
                        current_clust2 = 0
                if level2clust != current_clust2:
                    self.level2clustflow.append(0.)
                    current_clust2 += 1
                    self.tot_subclusters += 1
                    self.level2clusts[self.tot_subclusters] = []
                    # print "\nNew second cluster level:\n", self.tot_subclusters, level2clust, current_clust2
                    # print "Line:", line
                ##
                # print level1clust, level2clust
                # print current_clust1, current_clust2
                ##
                internal_flow = float(line[1])
                orig_node_name = int(line[3])
                self.weights.append(internal_flow)
                if self.read_clust_level==1: self.level1clusts[current_clust1].append(orig_node_name)
                else: self.level1clusts[self.tot_mainclusters].append(orig_node_name)
                self.level2clusts[self.tot_subclusters].append(orig_node_name)
                if self.read_clust_level==1: self.level1clustflow[current_clust1-1] += internal_flow
                else: self.level1clustflow[self.tot_mainclusters-1] += internal_flow
                self.level2clustflow[self.tot_subclusters-1] += internal_flow
                self.node_clust1_id[orig_node_name-1] = current_clust1
                self.node_clust2_id[orig_node_name-1] = self.tot_subclusters
                cum_flow += internal_flow
        ##quit("I quit")

        ##print "cumulative flow:", cum_flow
        assert abs(cum_flow-1.) < 1.0E-06, "error: tot. flow is not equal to unity"
        assert abs(sum(self.level1clustflow)-1.) < 1.0E-06, "error: tot. flow in first level partition" + \
                    " is not equal to unity"
        assert abs(sum(self.level2clustflow)-1.) < 1.0E-06, "error: tot. flow in second level" + \
                    " partition is not equal to unity"

        print "Number of level 1 clusters:", len(self.level1clustflow)
        print "Number of nodes in each level 1 cluster:\n", \
              [len(item) for key, item in self.level1clusts.iteritems()]
        print "Internal flow within level 1 clusters:\n", self.level1clustflow
        print "Total number of nodes in level 1 clusters:", \
              sum([len(item) for key, item in self.level1clusts.iteritems()])
        print "\n\nNumber of level 2 clusters:", len(self.level2clustflow)
        print "Number of nodes in each level 2 cluster:", \
              [len(item) for key, item in self.level2clusts.iteritems()]
        print "Internal flow within level 2 clusters:\n", self.level2clustflow
        print "Total number of nodes in level 2 clusters:", \
              sum([len(item) for key, item in self.level2clusts.iteritems()])

        #print "weights:\n", self.weights
        #print "level1clusts:\n", self.level1clusts
        #print "level2clusts:\n", self.level2clusts
        #print "level1clustflow:\n", self.level1clustflow
        #print "level2clustflow:\n", self.level2clustflow
        #print "node_clust1_id:\n", self.node_clust1_id
        #print "node_clust2_id:\n", self.node_clust2_id

        if not self.readflow: return

        # parse .flow file to calculate total flow between level 1 clusters or level 2 clusters
        # subcluster flow as adjacency matrix
        # NB the subcluster adjacency matrix has "self-loops" (internal flow within subcluster)
        if self.plot_level==1: # plot "level one" clusters
            n_subclusters = len(self.level1clustflow)
            node_clust_id = self.node_clust1_id
        elif self.plot_level==2: # plot "level two" clusters
            n_subclusters = self.tot_subclusters
            node_clust_id = self.node_clust2_id

        self.subclustflow = np.zeros((n_subclusters,n_subclusters),dtype=float)
        with open(self.rootfilename+".flow","r") as flowf:
            for line in flowf.readlines():
                line = line.split()
                if line[0] != "-->" and line[0] != "<--":
                    node1 = int(line[0])
                    level2clust1 = node_clust_id[node1-1]
                    continue
                elif line[0] == "<--": continue # we are just looking at outflows, not inflows
                node2 = int(line[1])
                flow_val = float(filter(lambda ch: ch not in "()", line[2]))
                level2clust2 = node_clust_id[node2-2]
                self.subclustflow[level2clust1-1,level2clust2-1] += flow_val

    ''' Function to plot the network and colour the communities '''
    def plot_network(self,clust_level,node_scale="uniform",label_nodes_bool=True, \
                     plot_major_flow_only=True):

        if clust_level==1: # colour by "level 1" clusters
            clust_id = self.node_clust1_id
            n_clusters = len(self.level1clusts)
        elif clust_level==2: # colour by "level 2" clusters
            clust_id = self.node_clust2_id
            n_clusters = self.tot_subclusters
        if not plot_major_flow_only or self.plot_level==0:
            colour_vals = [0. + float(i)*(1./float(n_clusters)) for i in range(n_clusters)]

        if self.plot_level==0: # plot individual nodes
            nodelist = [i+1 for i in range(self.n_nodes)]
            # set node colours according to clusters
            node_colours = [0. for i in range(self.n_nodes)]
            for i in range(self.n_nodes):
                node_colours[i] = colour_vals[clust_id[i]-1]

            # plot network
            vis_scale = 20000. # scale factor for size of nodes
            G = nx.DiGraph() # null directed graph
            G.add_nodes_from(nodelist) # add nodes to graph
            with open(self.rootfilename+".txt","r") as llf: # graph as linked list (Infomap input) file
                for line in llf.readlines():
                    line = line.split()
                    if line[0] == "#": continue # comment line
                    # add (weighted+directed) edges to graph sequentially
                    G.add_edge(int(line[0]),int(line[1]),weight=float(line[2]))
            plt.figure()
            nx.draw_networkx(G,with_labels=True,nodelist=nodelist,node_color=node_colours, \
                             node_size=[self.weights[i]*vis_scale for i in range(self.n_nodes)], \
                             cmap="hsv",vmin=0.,vmax=1.)
            plt.axis("off")
            plt.show()

        elif self.plot_level==1 or self.plot_level==2: # plot level 1 or 2 clusters (+ sub-cluster flow)

            if self.plot_level==1:
                n_subclusters = len(self.level1clustflow)
                levelclustflow = self.level1clustflow
            elif self.plot_level==2:
                n_subclusters = self.tot_subclusters
                levelclustflow = self.level2clustflow

            if plot_major_flow_only: # only plot subclusters with majority of the flow
                clust_id_sortbyflow_cut, clust_flow_vals_cut = self.find_major_flow()
                print "\nMajority of flow contained within first %i clusters" % len(clust_id_sortbyflow_cut)
                # nodelist = clust_id_sortbyflow_cut
                nodelist = [i for i in range(len(clust_id_sortbyflow_cut))]
                colour_vals = [0. + float(i)*(1./float(len(clust_id_sortbyflow_cut))) for i in \
                               range(len(clust_id_sortbyflow_cut))]
                self.clust_id_sortbyid = sorted(clust_id_sortbyflow_cut)
                node_colours = [colour_vals[self.clust_id_sortbyid.index(clust_id)] for clust_id in \
                                clust_id_sortbyflow_cut]

                # need to use a mask on self.subclustflow to make an adjacency matrix consisting only of
                # retained nodes (clusters with major flow)
                subclustflow_adjmtx = np.zeros((len(clust_id_sortbyflow_cut),len(clust_id_sortbyflow_cut)), \
                                               dtype=float)
                for i in range(len(clust_id_sortbyflow_cut)):
                    for j in range(len(clust_id_sortbyflow_cut)):
                        subclustflow_adjmtx[i,j] = self.subclustflow[clust_id_sortbyflow_cut[i]-1, \
                                                   clust_id_sortbyflow_cut[j]-1]
            else:
                nodelist = [i for i in range(n_subclusters)]
                node_colours = [colour_vals[i] for i in range(n_subclusters)]
                subclustflow_adjmtx = self.subclustflow

            vis_scale = 10000. # scale factor for size of nodes
            edge_vis_scale = 100. # scale factor for width of edges
            transparency_factor = 0.6 # edge transparency
            if node_scale == "uniform":
                node_size = 1.0E-02*vis_scale
            elif node_scale == "non-uniform":
                if plot_major_flow_only:
                    node_size=[clust_flow_vals_cut[i]*vis_scale for i in range(len(clust_id_sortbyflow_cut))]
                else:
                    node_size=[levelclustflow[i]*vis_scale for i in range(n_subclusters)]

            # plot network
            G = nx.DiGraph(subclustflow_adjmtx) # directed graph constructed from adjacency matrix as 2d np.array
            pos = nx.layout.spring_layout(G)
            plt.figure()
            nx.draw_networkx_nodes(G,pos,nodelist=nodelist,node_color=node_colours, \
                             node_size=node_size,cmap="hsv",vmin=0.,vmax=1.)
            edges = G.edges()
            if node_scale == "uniform":
                edge_weights = [1.0E-04]*len(edges)
            elif node_scale == "non-uniform":
                edge_weights = [0.]*len(edges)
                for i, edge in enumerate(edges):
                    edge_weights[i] = subclustflow_adjmtx[edge[0],edge[1]]
            nx.draw_networkx_edges(G,pos,width=[edge_weights[i]*edge_vis_scale for i in range(len(edges))], \
                                   alpha=transparency_factor,arrowstyle="->",arrowsize=20)
            if label_nodes_bool:
                #node_labels_dict = nx.get_node_attributes(G,"quack")
                #print "node_labels_dict:\n", node_labels_dict
                if plot_major_flow_only:
                    node_labels_dict = {i: clust_id for i, clust_id in enumerate(clust_id_sortbyflow_cut)}
                else:
                    node_labels_dict = {i-1: i for i in range(1,n_subclusters+1)}
                nx.draw_networkx_labels(G,pos,labels=node_labels_dict,font_size=10)
            plt.axis("off")
            plt.show()

    ''' Write group file (one entry per line, corresponding to node) e.g. for colouring disconn graph  '''
    def write_group_file(self):
        if self.plot_level==1: node_clust_id = self.node_clust1_id
        elif self.plot_level==2: node_clust_id = self.node_clust2_id
        new_id_vals = {clust_id: i+1 for i, clust_id in enumerate(self.clust_id_sortbyid)}
        ## TEST! (to make colouring easier to see) (HG-WC database) (second cluster level)
        ## new_id_vals = {1: 1, 2: 2, 3: 3, 4: 4, 8: 5, 9: 6, 10: 7, 11: 8}
        print "\nInfomap to disconn graph label mapping:\n", new_id_vals
        with open("groups.dat","w") as opf:
            for id_val in node_clust_id:
                if id_val not in self.clust_id_sortbyid:
                    opf.write("-1\n")
                else:
                    opf.write(str(new_id_vals[id_val]-1)+"\n")
                '''
                # TEST!
                if id_val not in [1, 2, 3, 4, 8, 9, 10, 11]:
                    opf.write("9\n")
                else:
                    opf.write(str(new_id_vals[id_val])+"\n")
                '''

    ''' Function to find subcluster IDs that give non-negligible (according to a flow_cutoff) contribution to
        the total flow. Use e.g. flow_cutoff=1.01 to return all clusters '''
    def find_major_flow(self,flow_cutoff=1.01):
        if self.plot_level==1: clust_flow_vals = self.level1clustflow
        elif self.plot_level==2: clust_flow_vals = self.level2clustflow
        clust_id_sortbyflow = sorted([i for i in range(1,len(clust_flow_vals)+1)],key = \
                                     lambda x: clust_flow_vals[x-1], reverse=True)
        # print clust_id_sortbyflow
        # print sorted(clust_flow_vals,reverse=True)
        cum_flow = 0. # cumulative flow
        for clust_id_cutoff, clust_id in enumerate(clust_id_sortbyflow):
            cum_flow += clust_flow_vals[clust_id-1]
            if cum_flow >= flow_cutoff: return clust_id_sortbyflow[:clust_id_cutoff], \
                                               sorted(clust_flow_vals,reverse=True)[:clust_id_cutoff]
        return clust_id_sortbyflow, sorted(clust_flow_vals,reverse=True)

if __name__=="__main__":
    import sys

    filename = sys.argv[1] # root file name for .tree and .txt files
    n_nodes = int(sys.argv[2]) # no. of nodes (see top of .log file)
    treeformat = int(sys.argv[3]) # read .tree (1) or .ftree (2) format
    readflow_int = int(sys.argv[4]) # read a .flow file (1/2 for Y/N)
    read_clust_level = int(sys.argv[5]) # (>=1). Read this cluster level + next cluster level
    clust_level = int(sys.argv[6]) # colour by cluster level (1,2)
    plot_level = int(sys.argv[7]) # plot indiv nodes (0) or plot clusters / "subclusters" (1/2)
    assert clust_level >= plot_level

    if readflow_int==1: readflow = True
    elif readflow_int==2: readflow = False

    parse_treefile1 = Parse_treefile(filename,n_nodes,treeformat,read_clust_level,readflow,plot_level)
            #pos = nx.get_node_attributes(G,"pos")
    parse_treefile1.plot_network(clust_level,node_scale="non-uniform")
    if parse_treefile1.plot_level!=0: parse_treefile1.write_group_file()
