#!/usr/bin/bash

#       Runs slice tree (several versions) for the synthetic graphs
#       @Arlei Silva

source default.sh

rm $results_time\_st.dat
rm $results_time\_stbs_fast.dat
rm $results_time\_stbs_slow.dat
rm $results_time\_al.dat
rm $results_time\_wvp.dat
rm $results_time\_wvb.dat

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
    prefix=$graph_name_prefix\_$g\_$p
    postfix=$g\_$p
    ctime=`grep compression_time out_st_$postfix.txt | cut -d ' ' -f3`
    avg_st=`echo "scale=10; $avg_st+$ctime" | bc`
    
    ctime=`grep compression_time out_al_$postfix.txt | cut -d ' ' -f3`
    avg_al=`echo "scale=10; $avg_al+$ctime" | bc`
    
    ctime=`grep compression_time out_wvp_$postfix.txt | cut -d ' ' -f3`
    avg_wvp=`echo "scale=10; $avg_wvp+$ctime" | bc`
    
    ctime=`grep compression_time out_wvb_$postfix.txt | cut -d ' ' -f3`
    avg_wvb=`echo "scale=10; $avg_wvb+$ctime" | bc`
    
    for((r=1; r<=$num_runs_sampling; r++))
    do
        postfix_samp=$postfix\_$r
        ctime=`grep compression_time out_stbs_fast_$postfix_samp.txt | cut -d ' ' -f3`
        avg_stbs_fast=`echo "scale=10; $avg_stbs_fast+$ctime" | bc`
        
	ctime=`grep compression_time out_stbs_slow_$postfix_samp.txt | cut -d ' ' -f3`
        avg_stbs_slow=`echo "scale=10; $avg_stbs_slow+$ctime" | bc`
    done
  done
  avg_stbs_fast=`echo "scale=10; $avg_stbs_fast/($num_graphs*$num_runs_sampling)" | bc`
  avg_stbs_slow=`echo "scale=10; $avg_stbs_slow/($num_graphs*$num_runs_sampling)" | bc`

  echo "$p	$avg_stbs_fast" >> $results_time\_stbs_fast.dat
  echo "$p	$avg_stbs_slow" >> $results_time\_stbs_slow.dat
  
  avg_st=`echo "scale=10; $avg_st/($num_graphs)" | bc`
  avg_al=`echo "scale=10; $avg_al/($num_graphs)" | bc`
  avg_wvp=`echo "scale=10; $avg_wvp/($num_graphs)" | bc`
  avg_wvb=`echo "scale=10; $avg_wvb/($num_graphs)" | bc`

  echo "$p	$avg_st" >> $results_time\_st.dat
  echo "$p	$avg_al" >> $results_time\_al.dat
  echo "$p	$avg_wvp" >> $results_time\_wvp.dat
  echo "$p	$avg_wvb" >> $results_time\_wvb.dat
done

