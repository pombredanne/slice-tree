#!/usr/bin/bash

#       Runs slice tree (several versions) for the synthetic graphs
#       @Arlei Silva

source default.sh

rm $results_sse_reduction.dat

for p in ${param_budget[@]}
do
  avg_wvp=0
  avg_wvb=0
  avg_stbs_fast=0
  avg_stbs_slow=0
  
  for((g=1; g<=$num_graphs; g++))
  do
    prefix=$graph_name_prefix\_$g\_$p
    postfix=$g\_$p
        
    alg_reduction=`grep sse_reduction out_wvp_$postfix.txt | cut -d ' ' -f3`
    alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
    sse=`grep -m1 sse out_wvp_$postfix.txt | cut -d ' ' -f3`
    alg_reduction=`echo "scale=10; $alg_reduction/($alg_reduction+$sse)" | bc`
	
    avg_wvp=`echo "scale=10; $avg_wvp+$alg_reduction" | bc`
    
    alg_reduction=`grep sse_reduction out_wvb_$postfix.txt | cut -d ' ' -f3`
    alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
    sse=`grep -m1 sse out_wvb_$postfix.txt | cut -d ' ' -f3`
    alg_reduction=`echo "scale=10; $alg_reduction/($alg_reduction+$sse)" | bc`

    avg_wvb=`echo "scale=10; $avg_wvb+$alg_reduction" | bc`

    for((r=1; r<=$num_runs_sampling; r++))
    do
        postfix_samp=$postfix\_$r
        alg_reduction=`grep sse_reduction out_stbs_fast_$postfix_samp.txt | cut -d ' ' -f3`
        alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
    	sse=`grep -m1 sse out_stbs_fast_$postfix_samp.txt | cut -d ' ' -f3`
    	alg_reduction=`echo "scale=10; $alg_reduction/($alg_reduction+$sse)" | bc`

	avg_stbs_fast=`echo "scale=10; $avg_stbs_fast+$alg_reduction" | bc`

        alg_reduction=`grep sse_reduction out_stbs_slow_$postfix_samp.txt | cut -d ' ' -f3`
        alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
    	sse=`grep -m1 sse out_stbs_slow_$postfix_samp.txt | cut -d ' ' -f3`
    	alg_reduction=`echo "scale=10; $alg_reduction/($alg_reduction+$sse)" | bc`
        
        avg_stbs_slow=`echo "scale=10; $avg_stbs_slow+$alg_reduction" | bc`
    done
  done
  
  avg_stbs_fast=`echo "scale=10; $avg_stbs_fast/($num_graphs*$num_runs_sampling)" | bc`
  avg_stbs_slow=`echo "scale=10; $avg_stbs_slow/($num_graphs*$num_runs_sampling)" | bc`
  avg_wvp=`echo "scale=10; $avg_wvp/($num_graphs)" | bc`
  avg_wvb=`echo "scale=10; $avg_wvb/($num_graphs)" | bc`

  echo "$p      $avg_wvb        $avg_wvp        $avg_stbs_slow  $avg_stbs_fast" >> $results_sse_reduction.dat
done
  

