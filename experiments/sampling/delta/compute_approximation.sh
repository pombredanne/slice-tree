#!/usr/bin/bash

#       Runs slice tree (several versions) for the synthetic graphs
#       @Arlei Silva

source default.sh

rm $results_approximation\_st.dat
rm $results_approximation\_stbs_fast.dat
rm $results_approximation\_stbs_slow.dat
rm $results_approximation\_stus_slow.dat

avg_st=0

for p in ${param_delta[@]}
do
  avg_stbs_fast=0
  avg_stbs_slow=0
  avg_stus_slow=0
  
  for((g=1; g<=$num_graphs; g++))
  do
    prefix=$graph_name_prefix\_$g
    postfix=$g
    optimal_reduction=`grep optimal_sse_reduction $prefix.stats | cut -d ' ' -f3`
    optimal_reduction=`echo ${optimal_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
    alg_reduction=`grep sse_reduction out_st_$postfix.txt | cut -d ' ' -f3`
    alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
    approximation=`echo "scale=10; $alg_reduction/$optimal_reduction" | bc`

    if [ "$approximation" > "1" ]
    then
      approximation=1
    fi

    avg_st=`echo "scale=10; $avg_st+" $approximation| bc`
    
    for((r=1; r<=$num_runs_sampling; r++))
    do
        postfix_samp=$postfix\_$p\_$r
        alg_reduction=`grep sse_reduction out_stbs_fast_$postfix_samp.txt | cut -d ' ' -f3`
        alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
        approximation=`echo "scale=10; $alg_reduction/$optimal_reduction" | bc`
	
        if [ "$approximation" > "1" ]
        then
          approximation=1
        fi

	avg_stbs_fast=`echo "scale=10; $avg_stbs_fast+$approximation" | bc`

        alg_reduction=`grep sse_reduction out_stbs_slow_$postfix_samp.txt | cut -d ' ' -f3`
        alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
        approximation=`echo "scale=10; $alg_reduction/$optimal_reduction" | bc`
	
        if [ "$approximation" > "1" ]
        then
          approximation=1
        fi
        
        avg_stbs_slow=`echo "scale=10; $avg_stbs_slow+$approximation" | bc`
        
        alg_reduction=`grep sse_reduction out_stus_slow_$postfix_samp.txt | cut -d ' ' -f3`
        alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
        approximation=`echo "scale=10; $alg_reduction/$optimal_reduction" | bc`
	
        if [ "$approximation" > "1" ]
        then
          approximation=1
	fi
	
        avg_stus_slow=`echo "scale=10; $avg_stus_slow+$approximation" | bc`
    done
  done
  avg_stbs_fast=`echo "scale=10; $avg_stbs_fast/($num_graphs*$num_runs_sampling)" | bc`
  avg_stbs_slow=`echo "scale=10; $avg_stbs_slow/($num_graphs*$num_runs_sampling)" | bc`
  avg_stus_slow=`echo "scale=10; $avg_stus_slow/($num_graphs*$num_runs_sampling)" | bc`

  echo "$p	$avg_stbs_fast" >> $results_approximation\_stbs_fast.dat
  echo "$p	$avg_stbs_slow" >> $results_approximation\_stbs_slow.dat
  echo "$p	$avg_stus_slow" >> $results_approximation\_stus_slow.dat
done
  
avg_st=`echo "scale=10; $avg_st/($num_graphs*${#param_delta[@]})" | bc`

for p in ${param_delta[@]}
do
  echo "$p	$avg_st" >> $results_approximation\_st.dat
done

