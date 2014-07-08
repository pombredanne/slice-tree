#!/bin/bash

#       Runs slice tree (several versions) for the synthetic graphs
#       @Arlei Silva

source default.sh

rm $results_sse_reduction.dat

echo "data      ST      STBS_SLOW       STBS_FAST" >> $results_sse_reduction.dat

for f in ${dblp_data_files[@]}
do
  avg_st=0
  avg_stbs_fast=0
  avg_stbs_slow=0
  avg_stus_slow=0
  postfix=$f  
  sse=`grep -m 1 sse out_st_$postfix.txt | cut -d ' ' -f3`
  sse=`echo ${sse} | sed -e 's/[eE]+*/\\*10\\^/'`
  alg_reduction=`grep sse_reduction out_st_$postfix.txt | cut -d ' ' -f3`
  alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`

  avg_st=`echo "scale=10; $alg_reduction"| bc`
  
  for((r=1; r<=$num_runs_sampling; r++))
  do
        postfix_samp=$postfix\_$r
        alg_reduction=`grep sse_reduction out_stbs_fast_$postfix_samp.txt | cut -d ' ' -f3`
        alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`

	avg_stbs_fast=`echo "scale=10; $avg_stbs_fast+$alg_reduction" | bc`

        alg_reduction=`grep sse_reduction out_stbs_slow_$postfix_samp.txt | cut -d ' ' -f3`
        alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
        
        avg_stbs_slow=`echo "scale=10; $avg_stbs_slow+$alg_reduction" | bc`
  done
  avg_stbs_fast=`echo "scale=10; $avg_stbs_fast/$num_runs_sampling" | bc`
  avg_stbs_slow=`echo "scale=10; $avg_stbs_slow/$num_runs_sampling" | bc`
  rate_st=`echo "scale=10; $avg_st/($sse)" | bc`
  rate_stbs_slow=`echo "scale=10; $avg_stbs_slow/($sse)" | bc`
  rate_stbs_fast=`echo "scale=10; $avg_stbs_fast/($sse)" | bc`
  echo "$f	$avg_st/$rate_st	$avg_stbs_slow/$rate_stbs_slow	$avg_stbs_fast/$rate_stbs_fast" >> $results_sse_reduction.dat
done
  


