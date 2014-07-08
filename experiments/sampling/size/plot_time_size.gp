set term postscript eps 40 enhanced color solid
set style fill solid 1.00 border
set encoding iso_8859_1
set output "time_size_syn.eps"
set xlabel "graph size"
set ylabel "time (sec.)"
set key top
set xrange [6000:2000000]
set yrange [:12000]
#set xtics 1
#set ytics 1000
set key top right
set format x "10^{%L}"
set format y "10^{%L}"
set log x
set log y
plot "size_time_st.dat" using 1:2 notitle with linespoints lt 1 lc 0 lw 4 pt 2 ps 4, "size_time_stus_slow.dat" using 1:2 notitle with linespoints lt 1 lc 2 lw 4 pt 5 ps 4,"size_time_stbs_slow.dat" using 1:2 notitle with linespoints lt 1 lc 3 lw 4 pt 7 ps 4,"size_time_stbs_fast.dat" using 1:2 notitle with linespoints lt 1 lc 1 lw 4 pt 9 ps 4

