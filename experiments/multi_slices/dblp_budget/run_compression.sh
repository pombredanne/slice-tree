#!/bin/bash

#       Runs slice tree (several versions) for the synthetic graphs
#       @Arlei Silva

source default.sh

for p in ${param_budget[@]}
do
    f=$data_file
    echo "$graph_compression -g dblp.graph -v $f.data -o output -c WVB -b $p > out_wvb_$p.txt"
    $graph_compression -g dblp.graph -v $f.data -o output -c WVB -b $p > out_wvb_$p.txt
    
    for((r=1; r<=$num_runs_sampling; r++))
    do
        postfix_samp=$p\_$r
        echo "$graph_compression -g dblp.graph -v $f.data -o output -c STBS -b $p -s dblp.sizes -m $max_radius -n $rate_fast_sampling -d $delta_fast_sampling -r $rho_fast_sampling -x > out_stbs_fast_$postfix_samp.txt"
        $graph_compression -g dblp.graph -v $f.data -o output -c STBS -b $p -s dblp.sizes -m $max_radius -n $rate_fast_sampling -d $delta_fast_sampling -r $rho_fast_sampling -x > out_stbs_fast_$postfix_samp.txt
        
	echo "$graph_compression -g dblp.graph -v $f.data -o output -c STBS -b $p -s dblp.sizes -m $max_radius -n $rate_slow_sampling -d $delta_slow_sampling -r $rho_slow_sampling -x > out_stbs_slow_$postfix_samp.txt"
        $graph_compression -g dblp.graph -v $f.data -o output -c STBS -b $p -s dblp.sizes -m $max_radius -n $rate_slow_sampling -d $delta_slow_sampling -r $rho_slow_sampling -x > out_stbs_slow_$postfix_samp.txt
    done
done
