import sys
import operator

num_snps = 8640
budget = 4923
relative_reductions = {}

for snp in range(0, num_snps):
  if snp % 12 == 0:
    file_name="out_slice_tree_"+str(snp)+"_"+str(budget)+".dat"
    snp_file = open(file_name, 'r')
    sse = 0
    sse_red = 0
    
    for line in snp_file:
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

for snp in sorted_reductions:
  file_name = snp[0]
  relative_reduction = snp[1]
  print "%s	%lf" % (file_name, relative_reduction) 

