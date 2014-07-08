#!/usr/bin/bash

#       Prepares ST and WV timing results for multi-slices for plotting
#       @Arlei Silva

source default.sh

for f in ${dblp_data_files[@]}
do
  rm ${f}_$results_time.dat
  echo "slices	WVP	ST	STI	STIF" >> ${f}_$results_time.dat
  for p in ${param_num_partitions[@]}
  do
    avg_wvp=0
    avg_st=0
    avg_stbs_fast=0
    avg_stbs_slow=0
    postfix=${f}_${p}
    
    ctime=`grep compression_time out_st_$postfix | cut -d ' ' -f3`
    avg_st=`echo "scale=10; $avg_st+$ctime" | bc`
    if [ -f out_wvp_$postfix ];
      then
      ctime=`grep compression_time out_wvp_$postfix | cut -d ' ' -f3`
      avg_wvp=`echo "scale=10; $avg_wvp+$ctime" | bc`
    fi

    for((r=1; r<=$num_runs_sampling; r++))
    do
        postfix_samp=$postfix\_$r
	if [ -f out_stbs_fast_$postfix_samp ];
          then
          ctime=`grep compression_time out_stbs_fast_$postfix_samp | cut -d ' ' -f3`
          avg_stbs_fast=`echo "scale=10; $avg_stbs_fast+$ctime" | bc`
        fi
    
        if [ -f out_stbs_slow_$postfix_samp ];
          then
	  ctime=`grep compression_time out_stbs_slow_$postfix_samp | cut -d ' ' -f3`
          avg_stbs_slow=`echo "scale=10; $avg_stbs_slow+$ctime" | bc`
        fi
    done
    avg_wvp=`echo "scale=10; $avg_wvp" | bc`
    avg_st=`echo "scale=10; $avg_st" | bc`
    avg_stbs_fast=`echo "scale=10; $avg_stbs_fast/$num_runs_sampling" | bc`
    avg_stbs_slow=`echo "scale=10; $avg_stbs_slow/$num_runs_sampling" | bc`

    echo "$p	$avg_wvp	$avg_st	$avg_stbs_slow	$avg_stbs_fast" >> ${f}_$results_time.dat
  done
done
