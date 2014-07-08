#!/usr/bin/bash
#bash run_compression.sh
./compute_time.sh
./compute_sse_reduction.sh
./compute_pruning.sh

# plot
tail -5 dblp_dm_dblp_time.dat > dblp_dm_time.dat
gnuplot plot_time_num_slices.gp
epstopdf dblp_dm_time.eps
scp dblp_dm_time.eps dblp_dm_time.pdf ~/Dropbox/slice_tree/paper/fig/
tail -5 dblp_dm_dblp_sse_reduction.dat > dblp_dm_sse_red.dat
gnuplot plot_time_sse_reduction.gp
epstopdf dblp_dm_sse_red.eps
scp dblp_dm_sse_red.eps dblp_dm_sse_red.pdf ~/Dropbox/slice_tree/paper/fig/
