#!/usr/bin/bash

#       Generates indexes for the synthetic graphs
#       @Arlei Silva

source default.sh

echo "$graph_compression -g dblp.graph -v dblp_dm.data -s dblp.sizes -t $num_threads -m $max_radius"
$graph_compression -g dblp.graph -v dblp_dm.data -s dblp.sizes -t $num_threads -m $max_radius
