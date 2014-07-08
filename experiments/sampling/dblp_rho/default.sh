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
param_rho=(0 0.2 0.4 0.6 0.8 1)
dblp_data_files=('dblp_alg')
results_time='dblp_rho_time'
results_pruning='dblp_rho_pruning'
results_approximation='dblp_rho_approximation'

