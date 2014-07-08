#!/bin/bash

#       Generates indexes for the synthetic graphs
#       @Arlei Silva

source default.sh

for d in ${data_files[@]}
  echo "$graph_compression -g $d.graph -v $d.data -s $d.sizes -t 1"
  $graph_compression -g $d.graph -v $d.data -s $d.sizes -t 1
