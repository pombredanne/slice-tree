#!/usr/bin/bash

#       Runs slice tree (several versions) for the synthetic graphs
#       @Arlei Silva

source default.sh

rm $results_sse_reduction.dat

for f in ${data_files[@]}
do
  postfix=$f  
  st=`grep sse_reduction out_st_$postfix.txt | cut -d ' ' -f3`
  st=`echo ${st} | sed -e 's/[eE]+*/\\*10\\^/'`
  sse=`grep -m1 sse out_st_$postfix.txt | cut -d ' ' -f3`
  sse=`echo ${sse} | sed -e 's/[eE]+*/\\*10\\^/'`
  st=`echo "scale=10; $st/($sse+$st)" | bc`
  
  wvb=`grep sse_reduction out_wvb_$postfix.txt | cut -d ' ' -f3`
  wvb=`echo ${wvb} | sed -e 's/[eE]+*/\\*10\\^/'`
  sse=`grep -m1 sse out_wvb_$postfix.txt | cut -d ' ' -f3`
  sse=`echo ${sse} | sed -e 's/[eE]+*/\\*10\\^/'`
  wvb=`echo "scale=10; $wvb/($sse+$wvb)" | bc`
  
  al=`grep sse_reduction out_al_$postfix.txt | cut -d ' ' -f3`
  al=`echo ${al} | sed -e 's/[eE]+*/\\*10\\^/'`
  sse=`grep -m1 sse out_al_$postfix.txt | cut -d ' ' -f3`
  sse=`echo ${sse} | sed -e 's/[eE]+*/\\*10\\^/'`
  al=`echo "scale=10; $al/($sse+$al)" | bc`

  echo "$f	$wvb	$al	$st" >> $results_sse_reduction.dat
done
  
