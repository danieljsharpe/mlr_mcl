#!/bin/bash

# script to repeatedly apply MLR-MCL to a KTN and output the results
# make sure that the arguments to the mlr_mcl exec are set correctly
# usage:
# ./drive_mcl.sh <n_it>

n_it=$1
seed=1 # first seed
# arguments to mlr_mcl that come before and after the seed, and extra args
# for the post-processing run to calculate quality metrics
mlr_mcl_args1='48800 62414 1.15 0.5 15 3 10 1.E-06 1.E-14'
mlr_mcl_args2='800'
mlr_mcl_args_xtra='20 0 1'


# cleanup
rm communities.*dat attractors.*dat mlr_mcl.*.out
rm modularity.dat avgncut.dat conductance.dat sce_node.*dat sce_edge.*dat
i=0
while [ $i -lt $n_it ]; do
    echo "MLR-MCL iteration:   $((i+1))"
    ./mlr_mcl $mlr_mcl_args1 $seed $mlr_mcl_args2 > mlr_mcl.$i.out
    ./mlr_mcl $mlr_mcl_args1 $seed $mlr_mcl_args2 $mlr_mcl_args_xtra > mlr_mcl.quality.$i.out
    mv communities_new.dat communities.$i.dat
    mv attractors_new.dat attractors.$i.dat
    rm communities.dat attractors.dat
    mv sce_node.dat sce_node.$i.dat
    mv sce_edge.dat sce_edge.$i.dat
    let i=i+1
    ((seed++))
done
