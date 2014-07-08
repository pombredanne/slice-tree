#!/usr/bin/bash

#       Runs slice tree (several versions) for the synthetic graphs
#       @Arlei Silva

source default.sh

rm $results_time.dat

for p in ${param_budget[@]}
do
  avg_stbs_fast=0
  avg_stbs_slow=0
  postfix=$p
  ctime=`grep compression_time out_wvb_$postfix.txt | cut -d ' ' -f3`
  avg_wvb=`echo "scale=10; $ctime" | bc`
  
  for((r=1; r<=$num_runs_sampling; r++))
  do
        postfix_samp=$postfix\_$r
        ctime=`grep compression_time out_stbs_fast_$postfix_samp.txt | cut -d ' ' -f3`
        avg_stbs_fast=`echo "scale=10; $avg_stbs_fast+$ctime" | bc`
        
	ctime=`grep compression_time out_stbs_slow_$postfix_samp.txt | cut -d ' ' -f3`
        avg_stbs_slow=`echo "scale=10; $avg_stbs_slow+$ctime" | bc`
  done
  avg_stbs_fast=`echo "scale=10; $avg_stbs_fast/$num_runs_sampling" | bc`
  avg_stbs_slow=`echo "scale=10; $avg_stbs_slow/$num_runs_sampling" | bc`

  echo "$p	$avg_wvb	$avg_stbs_slow	$avg_stbs_fast" >> $results_time.dat
done

