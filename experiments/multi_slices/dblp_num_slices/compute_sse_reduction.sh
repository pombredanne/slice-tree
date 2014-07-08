#!/bin/bash

#       Runs slice tree (several versions) for the synthetic graphs
#       @Arlei Silva

source default.sh

for f in ${dblp_data_files[@]}
do

  rm ${f}_$results_sse_reduction.dat
  echo "slices	WVP	ST	STI	STIF" >> ${f}_$results_sse_reduction.dat
  for p in ${param_num_partitions[@]}
  do
    avg_wvp=0
    avg_st=0
    avg_stbs_fast=0
    avg_stbs_slow=0

    postfix=${f}_${p}
    sse=`grep -m 1 sse out_st_$postfix | cut -d ' ' -f3`
    sse=`echo ${sse} | sed -e 's/[eE]+*/\\*10\\^/'`
    alg_reduction=`grep sse_reduction out_st_$postfix | cut -d ' ' -f3`
    alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
    avg_st=`echo "scale=10; $alg_reduction"| bc`  
    if [ -f out_wvp_$postfix ];
      then
      sse=`grep -m 1 sse out_wvp_$postfix | cut -d ' ' -f3`
      sse=`echo ${sse} | sed -e 's/[eE]+*/\\*10\\^/'`
      alg_reduction=`grep sse_reduction out_wvp_$postfix | cut -d ' ' -f3`
      alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
      avg_wvp=`echo "scale=10; $alg_reduction"| bc`
    fi
    for((r=1; r<=$num_runs_sampling; r++))
    do
        postfix_samp=$postfix\_$r
        if [ -f out_stbs_fast_$postfix_samp ];
          then
          alg_reduction=`grep sse_reduction out_stbs_fast_$postfix_samp | cut -d ' ' -f3`
          alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
          avg_stbs_fast=`echo "scale=10; $avg_stbs_fast+$alg_reduction" | bc`
        fi

        if [ -f out_stbs_slow_$postfix_samp ];
          then
          alg_reduction=`grep sse_reduction out_stbs_slow_$postfix_samp | cut -d ' ' -f3`
          alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
        
          avg_stbs_slow=`echo "scale=10; $avg_stbs_slow+$alg_reduction" | bc`
        fi
    done
    avg_stbs_fast=`echo "scale=10; $avg_stbs_fast/$num_runs_sampling" | bc`
    avg_stbs_slow=`echo "scale=10; $avg_stbs_slow/$num_runs_sampling" | bc`
    #rate_wvp=`echo "scale=10; $avg_wvp/($sse)" | bc`
    #rate_stbs_slow=`echo "scale=10; $avg_stbs_slow/($sse)" | bc`
    #rate_stbs_fast=`echo "scale=10; $avg_stbs_fast/($sse)" | bc`
    echo "$p	$avg_wvp	$avg_st	$avg_stbs_slow	$avg_stbs_fast" >> ${f}_$results_sse_reduction.dat
  done
done
  


