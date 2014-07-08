#!/usr/bin/bash

#       Runs slice tree (several versions) for the synthetic graphs
#       @Arlei Silva

source default.sh

rm $results_approximation.dat

for p in ${param_reduction[@]}
do
  avg_stbs_fast=0
  avg_stbs_slow=0
  avg_stus_slow=0
  
  for((g=1; g<=$num_graphs; g++))
  do
    prefix=$graph_name_prefix\_$g\_$p
    postfix=$g\_$p
    optimal_reduction=`grep sse_reduction out_st_$postfix.txt | cut -d ' ' -f3`
    optimal_reduction=`echo ${optimal_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`

    for((r=1; r<=$num_runs_sampling; r++))
    do
        postfix_samp=$postfix\_$r
        alg_reduction=`grep sse_reduction out_stbs_fast_$postfix_samp.txt | cut -d ' ' -f3`
        alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
        approximation=`echo "scale=10; $alg_reduction/$optimal_reduction" | bc`
	
        if [ $(echo " $approximation > 1" | bc) -eq 1 ]
        then
          approximation=1
        fi

	avg_stbs_fast=`echo "scale=10; $avg_stbs_fast+$approximation" | bc`

        alg_reduction=`grep sse_reduction out_stbs_slow_$postfix_samp.txt | cut -d ' ' -f3`
        alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
        approximation=`echo "scale=10; $alg_reduction/$optimal_reduction" | bc`
	
        if [ $(echo " $approximation > 1" | bc) -eq 1 ]
        then
          approximation=1
        fi
        
        avg_stbs_slow=`echo "scale=10; $avg_stbs_slow+$approximation" | bc`
        
        alg_reduction=`grep sse_reduction out_stus_slow_$postfix_samp.txt | cut -d ' ' -f3`
        alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
        approximation=`echo "scale=10; $alg_reduction/$optimal_reduction" | bc`
	
        if [ $(echo " $approximation > 1" | bc) -eq 1 ]
        then
          approximation=1
	fi
	
        avg_stus_slow=`echo "scale=10; $avg_stus_slow+$approximation" | bc`
    done
  done
  avg_stbs_fast=`echo "scale=10; $avg_stbs_fast/($num_graphs*$num_runs_sampling)" | bc`
  avg_stbs_slow=`echo "scale=10; $avg_stbs_slow/($num_graphs*$num_runs_sampling)" | bc`
  avg_stus_slow=`echo "scale=10; $avg_stus_slow/($num_graphs*$num_runs_sampling)" | bc`

  echo "$p      $avg_stus_slow  $avg_stbs_slow  $avg_stbs_fast" >> $results_approximation.dat
done
  


