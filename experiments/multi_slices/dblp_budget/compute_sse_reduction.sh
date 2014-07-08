#!/bin/bash

#       Runs slice tree (several versions) for the synthetic graphs
#       @Arlei Silva

source default.sh

rm $results_sse_reduction.dat

for p in ${param_budget[@]}
do
  avg_stbs_fast=0
  avg_stbs_slow=0
  postfix=$p  
  wvb=`grep sse_reduction out_wvb_$postfix.txt | cut -d ' ' -f3`
  wvb=`echo ${wvb} | sed -e 's/[eE]+*/\\*10\\^/'`
  sse=`grep -m1 sse out_wvb_$postfix.txt | cut -d ' ' -f3`
  sse=`echo ${sse} | sed -e 's/[eE]+*/\\*10\\^/'`
  wvb=`echo "scale=10; $wvb/($sse+$wvb)" | bc`

  for((r=1; r<=$num_runs_sampling; r++))
  do
        postfix_samp=$postfix\_$r
        alg_reduction=`grep sse_reduction out_stbs_fast_$postfix_samp.txt | cut -d ' ' -f3`
        alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
	sse=`grep -m1 sse out_stbs_fast_$postfix_samp.txt | cut -d ' ' -f3`
	sse=`echo ${sse} | sed -e 's/[eE]+*/\\*10\\^/'`
	alg_reduction=`echo "scale=10; $alg_reduction/($sse+$alg_reduction)" | bc`

	avg_stbs_fast=`echo "scale=10; $avg_stbs_fast+$alg_reduction" | bc`

        alg_reduction=`grep sse_reduction out_stbs_slow_$postfix_samp.txt | cut -d ' ' -f3`
        alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
	sse=`grep -m1 sse out_stbs_fast_$postfix_samp.txt | cut -d ' ' -f3`
	sse=`echo ${sse} | sed -e 's/[eE]+*/\\*10\\^/'`
	alg_reduction=`echo "scale=10; $alg_reduction/($sse+$alg_reduction)" | bc`
        
        avg_stbs_slow=`echo "scale=10; $avg_stbs_slow+$alg_reduction" | bc`
  done
  avg_stbs_fast=`echo "scale=10; $avg_stbs_fast/$num_runs_sampling" | bc`
  avg_stbs_slow=`echo "scale=10; $avg_stbs_slow/$num_runs_sampling" | bc`
  
  echo "$p	$wvb	$avg_stbs_slow	$avg_stbs_fast" >> $results_sse_reduction.dat
done
  
