set term postscript eps 40 enhanced color solid
set style fill solid 1.00 border
set encoding iso_8859_1
set output "time_radius_syn.eps"
set xlabel "radius"
set ylabel "time (sec.)"
set key top
set xrange [0.9:3.1]
set yrange [:899]
set xtics 1
set ytics 150
#set key top right
plot "radius_time_st.dat" using 1:2 title 'ST' with linespoints lt 1 lc 0 lw 4 pt 2 ps 4, "radius_time_stus_slow.dat" using 1:2 title 'STU' with linespoints lt 1 lc 2 lw 4 pt 5 ps 4,"radius_time_stbs_slow.dat" using 1:2 title 'STI' with linespoints lt 1 lc 3 lw 4 pt 7 ps 4,"radius_time_stbs_fast.dat" using 1:2 title 'STIF' with linespoints lt 1 lc 1 lw 4 pt 9 ps 4

