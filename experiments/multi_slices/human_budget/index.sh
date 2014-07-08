#!/usr/bin/bash

#       Generates indexes for the synthetic graphs
#       @Arlei Silva

source default.sh

echo "$graph_compression -g  $data_file.graph -v $data_file.data -s $data_file -t 1"
$graph_compression -g  $data_file.graph -v $data_file.data -s $data_file -t 1
