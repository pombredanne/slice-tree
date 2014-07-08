set term postscript eps 40 enhanced color solid
set style fill solid 1.00 border
set encoding iso_8859_1
set output "time_num_slices_syn.eps"
set xlabel "# partitions"
set ylabel "time (sec.)"
set key top
set xrange [1.5:40]
set yrange [:170]
#set xtics 1
set ytics 25
set key top right
set log x 2
plot "num_slices_time_st.dat" using 1:2 notitle with linespoints lt 1 lc 0 lw 4 pt 2 ps 4, "num_slices_time_stus_slow.dat" using 1:2 notitle with linespoints lt 1 lc 2 lw 4 pt 5 ps 4,"num_slices_time_stbs_slow.dat" using 1:2 notitle with linespoints lt 1 lc 3 lw 4 pt 7 ps 4,"num_slices_time_stbs_fast.dat" using 1:2 notitle with linespoints lt 1 lc 1 lw 4 pt 9 ps 4

