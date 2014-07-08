#!/bin/bash

#       Runs slice tree (several versions) for the synthetic graphs
#       @Arlei Silva

source default.sh

for p in ${param_num_partitions[@]}
do
  for g in ${dblp_data_files[@]}
  do
    #echo "$graph_compression -c ST -x -p $p -m $max_radius -g $data/dblp.graph -s $data/dblp.index -v $data/$g.data > out_st_${g}_${p}"
	#$graph_compression -c ST -x -p $p -m $max_radius -g $data/dblp.graph -s $data/dblp.index -v $data/$g.data > out_st_${g}_${p}

	#echo "$graph_compression -c AL -x -p $p -m $max_radius -g $data/dblp.graph -s $data/dblp.index -v $data/$g.data > out_al_${g}_${p}"
	#$graph_compression -c AL -x -p $p -m $max_radius -g $data/dblp.graph -s $data/dblp.index -v $data/$g.data > out_al_${g}_${p}

	#echo "$graph_compression -c WVP -x -p $p -m $max_radius -g $data/dblp.graph -s $data/dblp.index -v $data/$g.data > out_wvp_${g}_${p}"
	#$graph_compression -c WVP -x -p $p -m $max_radius -g $data/dblp.graph -s $data/dblp.index -v $data/$g.data > out_wvp_${g}_${p}

	echo "$graph_compression -c WVB -x -p $p -m $max_radius -g $data/dblp.graph -s $data/dblp.index -v $data/$g.data > out_wvb_${g}_${p}"
	$graph_compression -c WVB -x -p $p -m $max_radius -g $data/dblp.graph -s $data/dblp.index -v $data/$g.data > out_wvb_${g}_${p}
    
    #for((r=1; r<=$num_runs_sampling; r++))
    #do
		
		#echo "$graph_compression -c STBS -x -p $p -m $max_radius -g $data/dblp.graph -s $data/dblp.index -v $data/$g.data -n $rate_fast_sampling -d $delta_fast_sampling -r $rho_fast_sampling > out_stbs_${g}_${p}_${r}"
		#$graph_compression -c STBS -x -p $p -m $max_radius -g $data/dblp.graph -s $data/dblp.index -v $data/$g.data -n $rate_fast_sampling -d $delta_fast_sampling -r $rho_fast_sampling > out_stbs_${g}_${p}_${r}

	#echo "$graph_compression -c STBS -x -p $p -m $max_radius -g $data/dblp.graph -s $data/dblp.index -v $data/$g.data -n $rate_slow_sampling -d $delta_slow_sampling -r $rho_slow_sampling > out_stbs_slow_${g}_${p}_${r}"
	#$graph_compression -c STBS -x -p $p -m $max_radius -g $data/dblp.graph -s $data/dblp.index -v $data/$g.data -n $rate_slow_sampling -d $delta_slow_sampling -r $rho_slow_sampling > out_stbs_slow_${g}_${p}_${r}

    #done
  done
done
