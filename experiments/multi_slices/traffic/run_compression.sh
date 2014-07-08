#!/usr/bin/bash

#       Runs slice tree (several versions) for the synthetic graphs
#       @Arlei Silva

source default.sh

for f in ${data_files[@]}
do
    echo "$graph_compression -g traffic.graph -v $f.data -o output -c WVB -b $budget_alg > out_wvb_$f.txt"
    $graph_compression -g traffic.graph -v $f.data -o output -c WVB -b $budget_alg > out_wvb_$f.txt
    
    echo "$graph_compression -g traffic.graph -v $f.data -o output -c ST -b $budget_alg > out_st_$f.txt"
    $graph_compression -g traffic.graph -v $f.data -o output -c ST -b $budget_alg -s traffic.sizes > out_st_$f.txt
    
    echo "$graph_compression -g traffic.graph -v $f.data -o output -c AL -b $budget_alg > out_al_$f.txt"
    $graph_compression -g traffic.graph -v $f.data -o output -c AL -b $budget_alg > out_al_$f.txt
done

