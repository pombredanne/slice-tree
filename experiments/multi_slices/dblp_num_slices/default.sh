#!/bin/bash

#       Default parameters for experiments
#       @Arlei Silva

num_runs_sampling=10
max_radius=2
graph_compression='../../../graph_compression'
data='/home/petko/data/experiment/DBLP'
#dblp_data_files=('dblp_ai' 'dblp_alg' 'dblp_arch' 'dblp_bio' 'dblp_dm' 'dblp_edu' 'dblp_graph' 'dblp_hci' 'dblp_net' 'dblp_os' 'dblp_para' 'dblp_pl' 'dblp_sec' 'dblp_se')
#dblp_data_files=('dblp_dm' 'dblp_net' 'dblp_pl')
dblp_data_files=('dblp_dm')
#param_num_partitions=(16 32)
param_num_partitions=(2 4 8 16 32)

delta_fast_sampling=0.4
delta_slow_sampling=0.1
rho_fast_sampling=0.1
rho_slow_sampling=0.9
rate_fast_sampling=0.001
rate_slow_sampling=0.02
num_threads=1

results_time="dblp_time"
results_pruning="dblp_pruning"
results_sse_reduction="dblp_sse_reduction"
