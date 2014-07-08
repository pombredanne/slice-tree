#!/usr/bin/bash

#       Runs slice tree (several versions) for the synthetic graphs
#       @Arlei Silva

source default.sh

for p in ${param_sizes[@]}
do
  for((g=1; g<=$num_graphs; g++))
  do
    prefix=$graph_name_prefix\_$g\_$p
    postfix=$g\_$p
    echo "$graph_compression -g $prefix.graph -v $prefix.data -o output -c ST -p $num_partitions_alg -s $prefix.sizes -m $max_radius -x > out_st_$postfix.txt"
    $graph_compression -g $prefix.graph -v $prefix.data -o output -c ST -p $num_partitions_alg -s $prefix.sizes -m $max_radius -x > out_st_$postfix.txt
    
    echo "$graph_compression -g $prefix.graph -v $prefix.data -o output -c AL -p $num_partitions_alg -s $prefix.sizes -m $max_radius -x > out_al_$postfix.txt"
    $graph_compression -g $prefix.graph -v $prefix.data -o output -c AL -p $num_partitions_alg -s $prefix.sizes -m $max_radius -x > out_al_$postfix.txt
    
    echo "$graph_compression -g $prefix.graph -v $prefix.data -o output -c WVP -p $num_partitions_alg -s $prefix.sizes -m $max_radius -x > out_wvp_$postfix.txt"
    $graph_compression -g $prefix.graph -v $prefix.data -o output -c WVP -p $num_partitions_alg -s $prefix.sizes -m $max_radius -x > out_wvp_$postfix.txt
    
    echo "$graph_compression -g $prefix.graph -v $prefix.data -o output -c WVB -p $num_partitions_alg -s $prefix.sizes -m $max_radius -x > out_wvb_$postfix.txt"
    $graph_compression -g $prefix.graph -v $prefix.data -o output -c WVB -p $num_partitions_alg -s $prefix.sizes -m $max_radius -x > out_wvb_$postfix.txt
    
    for((r=1; r<=$num_runs_sampling; r++))
    do
        postfix_samp=$postfix\_$r
        echo "$graph_compression -g $prefix.graph -v $prefix.data -o output -c STBS -p $num_partitions_alg -s $prefix.sizes -m $max_radius -n $rate_fast_sampling -d $delta_fast_sampling -r $rho_fast_sampling -x > out_stbs_fast_$postfix_samp.txt"
        $graph_compression -g $prefix.graph -v $prefix.data -o output -c STBS -p $num_partitions_alg -s $prefix.sizes -m $max_radius -n $rate_fast_sampling -d $delta_fast_sampling -r $rho_fast_sampling -x > out_stbs_fast_$postfix_samp.txt
        echo "$graph_compression -g $prefix.graph -v $prefix.data -o output -c STBS -p $num_partitions_alg -s $prefix.sizes -m $max_radius -n $rate_slow_sampling -d $delta_slow_sampling -r $rho_slow_sampling -x > out_stbs_slow_$postfix_samp.txt"
        $graph_compression -g $prefix.graph -v $prefix.data -o output -c STBS -p $num_partitions_alg -s $prefix.sizes -m $max_radius -n $rate_slow_sampling -d $delta_slow_sampling -r $rho_slow_sampling -x > out_stbs_slow_$postfix_samp.txt
        
	echo "$graph_compression -g $prefix.graph -v $prefix.data -o output -c STUS -p $num_partitions_alg -s $prefix.sizes -m $max_radius -n $rate_fast_sampling -d $delta_fast_sampling -r $rho_fast_sampling -x > out_stus_fast_$postfix_samp.txt"
	$graph_compression -g $prefix.graph -v $prefix.data -o output -c STUS -p $num_partitions_alg -s $prefix.sizes -m $max_radius -n $rate_fast_sampling -d $delta_fast_sampling -r $rho_fast_sampling -x > out_stus_fast_$postfix_samp.txt
        echo "$graph_compression -g $prefix.graph -v $prefix.data -o output -c STUS -p $num_partitions_alg -s $prefix.sizes -m $max_radius -n $rate_slow_sampling -d $delta_slow_sampling -r $rho_slow_sampling -x > out_stus_slow_$postfix.txt"
        $graph_compression -g $prefix.graph -v $prefix.data -o output -c STUS -p $num_partitions_alg -s $prefix.sizes -m $max_radius -n $rate_slow_sampling -d $delta_slow_sampling -r $rho_slow_sampling -x > out_stus_slow_$postfix_samp.txt
    done
  done
done
