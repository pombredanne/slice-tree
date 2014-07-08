#!/usr/bin/bash

#       Generates indexes for the synthetic graphs
#       @Arlei Silva

source default.sh

echo "$graph_compression -g traffic.graph -v traffic.data -s traffic.sizes -t 1"
$graph_compression -g traffic.graph -v traffic_0.data -s traffic.sizes -t 1
