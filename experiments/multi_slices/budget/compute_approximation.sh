#!/usr/bin/bash

#       Runs slice tree (several versions) for the synthetic graphs
#       @Arlei Silva

source default.sh

rm $results_approximation\_st.dat
rm $results_approximation\_stbs_fast.dat
rm $results_approximation\_stbs_slow.dat
rm $results_approximation\_al.dat
rm $results_approximation\_wvp.dat
rm $results_approximation\_wvb.dat

avg_st=0

for p in ${param_budget[@]}
do
  avg_stbs_fast=0
  avg_stbs_slow=0
  avg_al=0
  avg_wvp=0
  avg_wvb=0
  
  for((g=1; g<=$num_graphs; g++))
  do
    prefix=$graph_name_prefix\_$g
    postfix=$g\_$p
    optimal_reduction=`grep optimal_sse_reduction $prefix.stats | cut -d ' ' -f3`
    optimal_reduction=`echo ${optimal_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
    
    alg_reduction=`grep sse_reduction out_st_$postfix.txt | cut -d ' ' -f3`
    alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
    approximation=`echo "scale=10; $alg_reduction/$optimal_reduction" | bc`

    if [ $(echo " $approximation > 1" | bc) -eq 1 ]
    then
      approximation=1
    fi

    avg_st=`echo "scale=10; $avg_st+ $approximation" | bc`
    
    alg_reduction=`grep sse_reduction out_al_$postfix.txt | cut -d ' ' -f3`
    alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
    approximation=`echo "scale=10; $alg_reduction/$optimal_reduction" | bc`

    if [ $(echo " $approximation > 1" | bc) -eq 1 ]
    then
      approximation=1
    fi

    avg_al=`echo "scale=10; $avg_al+$approximation" | bc`
    
    alg_reduction=`grep sse_reduction out_wvp_$postfix.txt | cut -d ' ' -f3`
    alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
    approximation=`echo "scale=10; $alg_reduction/$optimal_reduction" | bc`

    if [ $(echo " $approximation > 1" | bc) -eq 1 ]
    then
      approximation=1
    fi

    avg_wvp=`echo "scale=10; $avg_wvp+$approximation" | bc`
    
    alg_reduction=`grep sse_reduction out_wvb_$postfix.txt | cut -d ' ' -f3`
    alg_reduction=`echo ${alg_reduction} | sed -e 's/[eE]+*/\\*10\\^/'`
    approximation=`echo "scale=10; $alg_reduction/$optimal_reduction" | bc`

    if [ $(echo " $approximation > 1" | bc) -eq 1 ]
    then
      approximation=1
    fi

    avg_wvb=`echo "scale=10; $avg_wvb+$approximation" | bc`
    
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
    done
  done
  avg_stbs_fast=`echo "scale=10; $avg_stbs_fast/($num_graphs*$num_runs_sampling)" | bc`
  avg_stbs_slow=`echo "scale=10; $avg_stbs_slow/($num_graphs*$num_runs_sampling)" | bc`

  echo "$p	$avg_stbs_fast" >> $results_approximation\_stbs_fast.dat
  echo "$p	$avg_stbs_slow" >> $results_approximation\_stbs_slow.dat
  
  avg_st=`echo "scale=10; $avg_st/($num_graphs)" | bc`
  avg_al=`echo "scale=10; $avg_al/($num_graphs)" | bc`
  avg_wvp=`echo "scale=10; $avg_wvp/($num_graphs)" | bc`
  avg_wvb=`echo "scale=10; $avg_wvb/($num_graphs)" | bc`

  echo "$p	$avg_st" >> $results_approximation\_st.dat
  echo "$p	$avg_al" >> $results_approximation\_al.dat
  echo "$p	$avg_wvp" >> $results_approximation\_wvp.dat
  echo "$p	$avg_wvb" >> $results_approximation\_wvb.dat
done

