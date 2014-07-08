set term postscript eps 40 enhanced color solid
set style fill solid 1.00 border
set encoding iso_8859_1
set output "dblp_dm_time.eps"
set xlabel "# partitions"
set ylabel "time (sec.)"
set xrange [1.5:40]
set yrange [:13200]
#set xtics 1
set ytics 4000
set key center right
set log x 2
plot "dblp_dm_time.dat" using 1:2 title "WV" with linespoints lt 1 lc 5 lw 4 pt 8 ps 4, "dblp_dm_time.dat" using 1:3 title "ST" with linespoints lt 1 lc 0 lw 4 pt 2 ps 4,"dblp_dm_time.dat" using 1:4 title "STI" with linespoints lt 1 lc 3 lw 4 pt 7 ps 4, "dblp_dm_time.dat" using 1:5 title "STIF" with linespoints lt 1 lc 1 lw 4 pt 9 ps 4

