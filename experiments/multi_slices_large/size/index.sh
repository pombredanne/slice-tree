#!/usr/bin/bash

#       Generates indexes for the synthetic graphs
#       @Arlei Silva

source default.sh

for p in ${param_sizes[@]}
do
    for((g=1; g<=$num_graphs; g++))
    do
    	prefix=$graph_name_prefix\_$g\_$p
        echo "$graph_compression -g $prefix.graph -v $prefix.data -s $prefix.sizes -t $num_threads -m $max_radius"
        $graph_compression -g $prefix.graph -v $prefix.data -s $prefix.sizes -t $num_threads -m $max_radius
    done
done
