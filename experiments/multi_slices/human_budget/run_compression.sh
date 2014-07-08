#!/usr/bin/bash

#       Runs slice tree (several versions) for the synthetic graphs
#       @Arlei Silva

source default.sh

for p in ${param_budget[@]}
do
    postfix=$p
    f=$data_file
    echo "$graph_compression -g $f.graph -v $f.data -o output -c WVB -b $p > out_wvb_$postfix.txt"
    $graph_compression -g $f.graph -v $f.data -o output -c WVB -b $p > out_wvb_$postfix.txt
    
    echo "$graph_compression -g $f.graph -v $f.data -o output -c ST -b $p > out_st_$postfix.txt"
    $graph_compression -g $f.graph -v $f.data -o output -c ST -b $p -s $f.sizes > out_st_$postfix.txt
    
    echo "$graph_compression -g $f.graph -v $f.data -o output -c AL -b $p > out_al_$postfix.txt"
    $graph_compression -g $f.graph -v $f.data -o output -c AL -b $p > out_al_$postfix.txt
done

