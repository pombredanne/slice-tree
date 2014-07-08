#!/bin/bash

#       Default parameters for experiments
#       @Arlei Silva

num_graphs=10
num_runs_sampling=1
max_radius=2
num_edges=5
num_partitions_data=32
num_partitions_alg=32
num_vertices=100000
sse=200000
sse_reduction=100000
graph_generator='../../../graph_generator.py'
graph_compression='../../../graph_compression'
graph_name_prefix='syn'
param_num_partitions=(2 4 8 16 32)
delta_fast_sampling=0.4
delta_slow_sampling=0.1
rho_fast_sampling=0.1
rho_slow_sampling=0.6
rate_fast_sampling=0.01
rate_slow_sampling=0.1
num_threads=1
results_time='num_slices_time'
results_sse_reduction='num_slices_sse_reduction'
