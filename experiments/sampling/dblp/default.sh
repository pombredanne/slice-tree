#!/bin/bash

#       Default parameters for experiments
#       @Arlei Silva

num_runs_sampling=10
max_radius=2
num_partitions_alg=2
graph_compression='../../../graph_compression'
delta_fast_sampling=0.4
delta_slow_sampling=0.1
rho_fast_sampling=0.1
rho_slow_sampling=0.9
rate_fast_sampling=0.001
rate_slow_sampling=0.02
num_threads=8
#dblp_data_files=('dblp_ai' 'dblp_alg' 'dblp_arch' 'dblp_bio' 'dblp_dm' 'dblp_edu' 'dblp_graph' 'dblp_hci' 'dblp_net' 'dblp_os' 'dblp_para' 'dblp_pl' 'dblp_sec' 'dblp_se')
dblp_data_files=('dblp_alg' 'dblp_sec' 'dblp_dm' 'dblp_net')
results_time="dblp_time"
results_pruning="dblp_pruning"
results_approximation="dblp_approximation"
