#!/bin/bash

#       Runs slice tree (several versions) for the synthetic graphs
#       @Arlei Silva

source default.sh

rm $results_approximation.dat

for p in ${param_rho[@]}
do
  avg_stbs_fast=0
  avg_stbs_slow=0
  
  postfix=${dblp_data_files[0]}
  optimal_reduction=`grep sse_reduction out_st_$postfix.txt | cut -d ' ' -f3`
  optimal_reduction=`echo ${optimal_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
  
  for((r=1; r<=$num_runs_sampling; r++))
  do
      postfix_samp=$postfix\_$p\_$r

      alg_reduction=`grep sse_reduction out_stbs_fast_$postfix_samp.txt | cut -d ' ' -f3`
      alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
      approximation=`echo "scale=10; $alg_reduction/$optimal_reduction" | bc`
      
      avg_stbs_fast=`echo "scale=10; $avg_stbs_fast+$approximation" | bc`

      alg_reduction=`grep sse_reduction out_stbs_slow_$postfix_samp.txt | cut -d ' ' -f3`
      alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
      approximation=`echo "scale=10; $alg_reduction/$optimal_reduction" | bc`

      avg_stbs_slow=`echo "scale=10; $avg_stbs_slow+$approximation" | bc`
  done
  avg_stbs_fast=`echo "scale=10; $avg_stbs_fast/($num_runs_sampling)" | bc`
  avg_stbs_slow=`echo "scale=10; $avg_stbs_slow/($num_runs_sampling)" | bc`

  echo "$p	$avg_stbs_slow	$avg_stbs_fast" >> $results_approximation.dat
done
  
