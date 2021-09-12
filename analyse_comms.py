'''
Python script to analyse communities of a kinetic transition network:
clustering entropy, quality functions, community sizes, scaled link density of communities
'''

from math import log
from math import sqrt

### SET PARAMS ###
n_nodes=48800
n_edges=62414
n_its=100

# CLUSTERING ENTROPY

ic_node_probs = [0.]*n_nodes # probability of a node being inter-community
ic_edge_probs = [0.]*n_edges # probability of an edge being inter-community

for i in range(n_its):
    with open("sce_node."+str(i)+".dat","r") as sce_node_f:
        for j, line in enumerate(sce_node_f.readlines()):
            ic_node_probs[j] += float(line.split()[0])
    with open("sce_edge."+str(i)+".dat","r") as sce_edge_f:
        for j, line in enumerate(sce_edge_f.readlines()):
            if int(line.split()[0])==-1: continue # dead TS flag
            ic_edge_probs[j] += float(line.split()[0])
for i in range(n_nodes):
    ic_node_probs[i] *= 1./float(n_its)
for i in range(n_edges):
    ic_edge_probs[i] *= 1./float(n_its)

s_ce = 0. # clustering entropy
for j in range(n_edges):
    p_ij = ic_edge_probs[j]
    if p_ij==0. or p_ij==1.: continue
    else: s_ce += (p_ij*log(p_ij,2.)) + ((1.-p_ij)*log(1.-p_ij,2.))
s_ce *= -1./float(n_edges)

print "clustering entropy is: %1.6f" % s_ce

with open("ic_node_probs.dat","w") as icnp_f:
    for j in range(n_nodes):
        icnp_f.write("%1.6f\n" % ic_node_probs[j])
with open("ic_edge_probs.dat","w") as icep_f:
    for j in range(n_edges):
        icep_f.write("%1.6f\n" % ic_edge_probs[j])

bin_vals = [0.+(float(i)*0.05) for i in range(20+1)]
ic_node_probs_bins, ic_edge_probs_bins = [0]*20, [0]*20
curr_bin_idx = 0
for p_i in sorted(ic_node_probs):
    while p_i > bin_vals[curr_bin_idx+1]: curr_bin_idx += 1
    ic_node_probs_bins[curr_bin_idx] += 1
curr_bin_idx = 0
for p_ij in sorted(ic_edge_probs):
    while p_ij > bin_vals[curr_bin_idx+1]: curr_bin_idx += 1
    ic_edge_probs_bins[curr_bin_idx] += 1

hist_files = ["ic_node_probs_distrib", "ic_edge_probs_distrib"]
bin_arrs = [ic_node_probs_bins,ic_edge_probs_bins]
for i in range(2):
    hist_file = hist_files[i]
    bin_arr = bin_arrs[i]
    with open(hist_file+".dat","w") as hist_f:
        for bin_idx, bin_popn in enumerate(bin_arr):
            hist_f.write("%1.2f   %5i\n" % (bin_vals[bin_idx],bin_popn))

quit()

# AVG & STD DEV OF QUALITY FUNCTIONS

quality_files = ["modularity","avgncut","conductance"]
for quality_f in quality_files:
    qf_vals = [0.]*n_its;
    qf_avg, qf_sd = 0., 0.
    with open(quality_f+".dat","r") as qual_f:
        for i in enumerate(qual_f.readlines()):
            qf_vals[i] = float(line.split()[0])
            qf_avg += qf_vals[i]
    qf_avg *= 1./float(n_its)
    for i in range(n_its):
        qf_sd += (qf_vals[i]-qf_avg)**2
    qf_sd = sqrt((1./float(n_its-1))*qf_sd)
    print "%s:    avg:  %1.8f    sd:  %1.8f" % (quality_f,qf_avg,qf_sd)

# COMMUNITY SIZES DISTRIBUTION

# SCALED LINK DENSITY OF COMMUNITIES DISTRIBUTION
