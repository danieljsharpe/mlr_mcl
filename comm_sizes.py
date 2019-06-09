''' Write out sizes of each community to file "comm_sizes.dat"
    Set condition for writing a file "pick.dat" (to show only certain nodes with disconnectionDPS)
    Set condition for writing a file "colours.dat" (colouring of communities)
    Read/write attractors, 1/2 for Y/N '''

import sys

ncomm = int(sys.argv[1])
nmin = int(sys.argv[2])
do_att = int(sys.argv[3])

comm_sizes = [0]*ncomm
comm_ids = [0]*nmin

### SET CONDITIONS ### (x is a min_id)
#pick_condition = lambda x: (comm_ids[x]==173)
pick_condition = lambda x: (comm_sizes[comm_ids[x]]>=100)
colour_condition = lambda x: (comm_sizes[comm_ids[x]]>=100)

# find communities of minima and community sizes
with open("communities.dat","r") as comm_f:
    min_id=0
    for line in comm_f.readlines():
        comm_ids[min_id] = int(line)
        comm_sizes[comm_ids[min_id]] += 1
        min_id += 1

# find attractors
if do_att==1:
    is_attractor = {i: False for i in range(nmin)}
    found_attractor = {i: False for i in range(ncomm)}
    with open("attractors.dat","r") as att_f:
        for line in att_f.readlines():
            is_attractor[int(line)-1] = True

# write community sizes to file
with open("comm_sizes.dat","w") as commsz_f:
    for comm_sz in comm_sizes:
        commsz_f.write("%i\n" % comm_sz)


pick_f = open("pick.dat","w")
colours_f = open("colours.dat","w")
if do_att==1: newatt_f = open("attractors_pick.dat","w")

min_valid = {i: (True if colour_condition(i) else False) for i in range(nmin)}
new_comm_ids = {i: -1 for i in range(ncomm)} # ALL invalid communities are assigned a new ID of -1 by default
nseen = 1

# write "pick.dat", "attractors_pick.dat" and "colours.dat"
for i in range(nmin):
    if pick_condition(i):
        pick_f.write("%i\n" % (i+1))
    if do_att==1:
        if is_attractor[i] and not found_attractor[comm_ids[i]] and colour_condition(i):
            newatt_f.write("%i\n" % (i+1))
            found_attractor[comm_ids[i]] = True
    if min_valid[i] and new_comm_ids[comm_ids[i]] == -1: # min is valid and we have not seen its community yet
        new_comm_ids[comm_ids[i]] = nseen # community given a new id
        nseen += 1
    colours_f.write("%i\n" % new_comm_ids[comm_ids[i]])
print "there are", nseen-1, "communities that obey the colour_condition"

pick_f.close()
colours_f.close()
if do_att==1: newatt_f.close()
