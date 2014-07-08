#!/bin/bash

#       Runs slice tree (several versions) for the synthetic graphs
#       @Arlei Silva

source default.sh

rm $results_time.dat

avg_st=0

for p in ${param_rho[@]}
do
  avg_stbs_fast=0
  avg_stbs_slow=0
  
  prefix=$graph_name_prefix\_$g\_$p
  postfix=${dblp_data_files[0]}
  ctime=`grep compression_time out_st_$postfix.txt | cut -d ' ' -f3`
  avg_st=`echo "scale=10; $ctime" | bc`
    
  for((r=1; r<=$num_runs_sampling; r++))
  do
      postfix_samp=$postfix\_$p\_$r
      ctime=`grep compression_time out_stbs_fast_$postfix_samp.txt | cut -d ' ' -f3`
      avg_stbs_fast=`echo "scale=10; $avg_stbs_fast+$ctime" | bc`
        
      ctime=`grep compression_time out_stbs_slow_$postfix_samp.txt | cut -d ' ' -f3`
      avg_stbs_slow=`echo "scale=10; $avg_stbs_slow+$ctime" | bc`
  done
  avg_stbs_fast=`echo "scale=10; $avg_stbs_fast/($num_runs_sampling)" | bc`
  avg_stbs_slow=`echo "scale=10; $avg_stbs_slow/($num_runs_sampling)" | bc`

  echo "$p	$avg_st	$avg_stbs_slow	$avg_stbs_fast" >> $results_time.dat
done
  
