#!/usr/bin/python2
import sys
import operator

num_tissues = 114
num_partitions = 32
relative_reductions = {}

for t in range(0, num_tissues):
  if t % 12 == 0:
    file_name="out_st_"+str(t)+"_"+str(num_partitions)+".dat"
    t_file = open(file_name, 'r')
    sse = 0
    sse_red = 0
    
    for line in t_file:
      line = line.rstrip()
      vec = line.rsplit()
      stat = vec[0]
      value = float(vec[2])

      if stat == "sse":
        sse = value
      if stat == "sse_reduction":
        sse_red = value
    
    relative_reductions[file_name] = float(sse_red) / (sse + sse_red)

sorted_reductions = sorted(relative_reductions.iteritems(), key=operator.itemgetter(1), reverse=True)

for t in sorted_reductions:
  file_name = t[0]
  relative_reduction = t[1]
  print "%s	%lf" % (file_name, relative_reduction) 

