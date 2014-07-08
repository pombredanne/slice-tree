set term postscript eps 40 enhanced color solid
set style fill solid 1.00 border
set style data histogram
set style histogram cluster gap 1
set encoding iso_8859_1
set output "pruning_radius_syn.eps"
set xlabel "radius"
set ylabel "% pruning #samples"
set key top
set yrange [0:1]
set xrange [-1.5:4.2]
set boxwidth .9 absolute
set key at 4.5, 1
set key samplen 1
plot "radius_pruning.dat" using 2:xtic(1) title 'STI' lc 3 fs pattern 4,'' using 3 title 'STIF' lc 1 fs pattern 1

