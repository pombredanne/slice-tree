#!/bin/bash

#       Default parameters for experiments
#       @Arlei Silva

num_runs_sampling=10
max_radius=2
param_budget=(25 50 100 200 400)
graph_compression='../../../graph_compression'
delta_fast_sampling=0.4
delta_slow_sampling=0.1
rho_fast_sampling=0.1
rho_slow_sampling=0.9
rate_fast_sampling=0.001
rate_slow_sampling=0.02
num_threads=8
data_file='dblp_alg'
results_time="dblp_time_budget"
results_pruning="dblp_pruning_budget"
results_sse_reduction="dblp_sse_reduction_budget"
