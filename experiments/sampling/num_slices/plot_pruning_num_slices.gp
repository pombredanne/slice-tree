set term postscript eps 40 enhanced color solid
set style fill solid 1.00 border
set style data histogram
set style histogram cluster gap 1
set encoding iso_8859_1
set output "pruning_num_slices_syn.eps"
set xlabel "# partitions"
set ylabel "% pruning #samples"
set key top
set yrange [0:1]
#set xrange [-1:3]
set boxwidth .9 absolute
set xtics ("2" 0, "4" 1, "8" 2, "16" 3, "32" 4)
plot "num_slices_pruning.dat" using 2 notitle lc 3 fs pattern 4,'' using 3 notitle lc 1 fs pattern 1

