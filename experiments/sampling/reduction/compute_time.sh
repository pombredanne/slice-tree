#!/usr/bin/bash

#       Runs slice tree (several versions) for the synthetic graphs
#       @Arlei Silva

source default.sh

rm $results_time\_st.dat
rm $results_time\_stbs_fast.dat
rm $results_time\_stbs_slow.dat
rm $results_time\_stus_slow.dat

for p in ${param_reduction[@]}
do
  avg_st=0
  avg_stbs_fast=0
  avg_stbs_slow=0
  avg_stus_slow=0
  for((g=1; g<=$num_graphs; g++))
  do
    prefix=$graph_name_prefix\_$g\_$p
    postfix=$g\_$p
    ctime=`grep compression_time out_st_$postfix.txt | cut -d ' ' -f3`
    avg_st=`echo "scale=10; $avg_st+$ctime" | bc`
    for((r=1; r<=$num_runs_sampling; r++))
    do
        postfix_samp=$postfix\_$r
        ctime=`grep compression_time out_stbs_fast_$postfix_samp.txt | cut -d ' ' -f3`
        avg_stbs_fast=`echo "scale=10; $avg_stbs_fast+$ctime" | bc`
        
	ctime=`grep compression_time out_stbs_slow_$postfix_samp.txt | cut -d ' ' -f3`
        avg_stbs_slow=`echo "scale=10; $avg_stbs_slow+$ctime" | bc`
        
	ctime=`grep compression_time out_stus_slow_$postfix_samp.txt | cut -d ' ' -f3`
        avg_stus_slow=`echo "scale=10; $avg_stus_slow+$ctime" | bc`
    done
  done
  avg_st=`echo "scale=10; $avg_st/$num_graphs" | bc`
  avg_stbs_fast=`echo "scale=10; $avg_stbs_fast/($num_graphs*$num_runs_sampling)" | bc`
  avg_stbs_slow=`echo "scale=10; $avg_stbs_slow/($num_graphs*$num_runs_sampling)" | bc`
  avg_stus_slow=`echo "scale=10; $avg_stus_slow/($num_graphs*$num_runs_sampling)" | bc`

  echo "$p	$avg_st" >> $results_time\_st.dat
  echo "$p	$avg_stbs_fast" >> $results_time\_stbs_fast.dat
  echo "$p	$avg_stbs_slow" >> $results_time\_stbs_slow.dat
  echo "$p	$avg_stus_slow" >> $results_time\_stus_slow.dat
done
