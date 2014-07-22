/**
Copyright (c) 2013, Arlei Silva
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

@author: Arlei Silva (arleilps@gmail.com)
**/

/**
 *  Implementation of graph compression
**/

/*std includes*/
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <climits>
#include <limits>
#include <cfloat>
#include <stdint.h>

/*my includes*/
#include "graph_compression.h"

double GraphCompression::sse_value = 0;
double GraphCompression::sse_reduction_value = 0;
double GraphCompression::compression_rate_value = 0;
unsigned int GraphCompression::budget_value = 0;
double GraphCompression::maximum_pointwise_error_value = 0;
double GraphCompression::peak_signal_to_noise_ratio_value = 0;
double GraphCompression::root_mean_squared_error_value = 0;
double GraphCompression::normalized_sse_value = 0;
GraphCompressionAlgorithm* GraphCompression::compression_algorithm = NULL;
unsigned int SliceTreeSamp::num_pruned_bound_1 = 0;
unsigned int SliceTreeSamp::num_pruned_bound_2 = 0;
unsigned int SliceTreeSamp::num_pruned_bound_3 = 0;
unsigned int SliceTreeSamp::num_pruned = 0;
unsigned int SliceTreeSamp::total_slices = 0;

/**
 * Computes epsilon, as defined for hoeffding bounds.
 * @num_samples number of samples
 * @delta probability that the bounds don't hold
 * @return value of epsilon
 * @throws
**/
double compute_epsilon_num_samples(
    unsigned int num_samples, double delta)
{
    if(num_samples > 0)
    {
        return sqrt(-1 * (double) log(delta) / (2 * num_samples));
    }
    else
    {
        return  std::numeric_limits<float>::max();
    }
}

/**
 * Computes epsilon, as defined for hoeffding bounds.
 * @param theta range for values
 * @num_samples number of samples
 * @delta probability that the bounds don't hold
 * @return value of epsilon
 * @throws
**/
double compute_epsilon_mean(double theta,
    unsigned int num_samples, double delta)
{
    if(num_samples > 0)
    {
        return sqrt(-1 * (double)(pow(theta, 2) *
            log((double) delta / 2)) / (2 * num_samples));
    }
    else
    {
        return  std::numeric_limits<float>::max();
    }
}
/**
 * Compresses the graph using the compression algorithm
 * @param graph_to_compress graph to be compressed
 * @param budget size of the the compressed file in bytes
 * @param output_file_name output file name
 * @return
 * @throws
**/
void GraphCompression::compress(Graph& graph_to_compress,
    GraphCompressionAlgorithm& algorithm, const unsigned int budget,
    const std::string& output_file_name)
{
    compression_algorithm = &algorithm;
    budget_value = budget;
    compression_algorithm->compress(budget);
    compression_algorithm->write(output_file_name);
    compression_algorithm->set_compressed_values();
    compute_statistics(graph_to_compress);
}

/**
 * Computes the sse of the compression
 * @param
 * @return sse
 * @throws
**/
void GraphCompression::compute_sse(Graph& graph)
{
    for(unsigned int v = 0; v < graph.size(); v++)
    {
        sse_value = sse_value + pow(graph.orig_value(v) - compression_algorithm->value(v), 2);
    }
}

/**
 * Computes the sse reduction of the compression
 * reduction = (sse data - sse compression) / sse data
 * @param
 * @return
 * @throws
**/
void GraphCompression::compute_sse_reduction(Graph& graph)
{
    double initial_sse = 0;
    double average = 0;

    for(unsigned int v = 0; v < graph.size(); v++)
    {
        average = average + graph.orig_value(v);
    }

    average = (double) average / graph.size();

    for(unsigned int v = 0; v < graph.size(); v++)
    {
        initial_sse = initial_sse + pow(average-graph.orig_value(v), 2);
    }

    sse_reduction_value = (double) (initial_sse - sse_value);
}

/**
 * Computes the normalized sse of the compression
 * normalized sse = sse / sum values squared
 * @param
 * @return
 * @throws
**/
void GraphCompression::compute_normalized_sse(Graph& graph)
{
    double sum_values_squared = 0;

    for(unsigned int v = 0; v < graph.size(); v++)
    {
        sum_values_squared = sum_values_squared + pow(graph.orig_value(v), 2);
    }

    normalized_sse_value = (double) sse_value / sum_values_squared;
}

/**
 * Computes the compression rate
 * compression rate = size original data / size compressed data
 * @param
 * @return
 * @throws
**/
void GraphCompression::compute_compression_rate(Graph& graph)
{
    double budget_graph = graph.size() * SIZE_FLOAT_INT;

    compression_rate_value = (double) budget_graph / budget_value;
}

/**
 * Computes the root mean squared error of the compression
 * RMSE = sqrt(sse/number of vertices)
 * @param
 * @return
 * @throws
**/
void GraphCompression::compute_root_mean_squared_error(Graph& graph)
{
    root_mean_squared_error_value = sqrt((double) sse_value / graph.size());
}

/**
 * Computes the maximum pointwise error of the compression
 * MPE = max difference between value and recovered value in absolute value
 * @param
 * @return
 * @throws
**/
void GraphCompression::compute_maximum_pointwise_error(Graph& graph)
{
    double pointwise_error;
    for(unsigned int v = 0; v < graph.size(); v++)
    {
        pointwise_error = fabs(graph.orig_value(v) - compression_algorithm->value(v));

        if(pointwise_error > maximum_pointwise_error_value)
        {
            maximum_pointwise_error_value = pointwise_error;
        }
    }
}

/**
 * Computes the peak signal to noise ratio of the compression
 * PSNR = 20*log_10(max value / RMSE)
 * @param
 * @return
 * @throws
**/
void GraphCompression::compute_peak_signal_to_noise_ratio(Graph& graph)
{
    double max_value = 0;

    for(unsigned int v = 0; v < graph.size(); v++)
    {
        if(graph.orig_value(v) > max_value)
        {
            max_value = graph.orig_value(v);
        }
    }

    if(root_mean_squared_error_value > 0)
    {
        peak_signal_to_noise_ratio_value =
            (double) 20 * log10((double) max_value
            / root_mean_squared_error_value);
    }
}

/**
 * Computes several statistics (mostly error measures) for
 * the compression
 * @param graph graph
 * @return
 * @throws
**/
void GraphCompression::compute_statistics(Graph& graph)
{
    compute_sse(graph);
    compute_sse_reduction(graph);
//  compute_normalized_sse(graph);
    compute_compression_rate(graph);
//  compute_root_mean_squared_error(graph);
//  compute_maximum_pointwise_error(graph);
//  compute_peak_signal_to_noise_ratio(graph);
}

/**
 * Decompresses the graph using the compression algorithm
 * @param compressed_file_name compressed file
 * @return the recovered graph
 * @throws
**/
void GraphCompression::decompress(const std::string& compressed_file_name,
    GraphCompressionAlgorithm& algorithm, Graph& graph)
{
    algorithm.decompress();
}

/**
 * Constructor.
 * @param graph_to_compress input graph
 * @param budget available budget
 * @return
 * @throws
**/
GraphCompressionAlgorithm::GraphCompressionAlgorithm(Graph& graph_to_compress)
{
    graph = &graph_to_compress;
    values.reserve(graph->size());
}

/**
 * Constructor.
 * @param input_file_name file withe the compressed data
 * @param graph_to_decompress graph to be decompressed
 * @return
 * @throws
**/
GraphCompressionAlgorithm::GraphCompressionAlgorithm(
    const std::string& input_file_name, Graph& graph_to_decompress)
{
    graph = &graph_to_decompress;
    compressed_file_name = input_file_name;
    values.reserve(graph->size());
}

/**
 * Identifies the optimal cut (center/radius) for the
 * partition represented as a slice tree node
 * @param st_node slice tree node
 * @return
 * @throws
**/
void SliceTree::optimal_cut(st_node_t* st_node)
{
    unsigned int best_center = st_node->partition[0];
    unsigned int best_radius = 0;
    double min_error = std::numeric_limits<double>::max();
    std::pair<double , unsigned int> e_r;
    double error;
    unsigned int radius;

    /*Computes the best radius for each possible center in the
    * partition*/
    for(unsigned int c = 0; c < st_node->partition.size(); c++)
    {
        e_r = min_error_radius(st_node->partition[c],
            st_node->partition, st_node->diameter,
            st_node->in_partition, st_node->average);
        error = global_error - e_r.first;
        radius = e_r.second;

        if(error < min_error)
        {
            best_center = st_node->partition[c];
            best_radius = radius;
            min_error = error;
        }
    }

    st_node->center = best_center;
    st_node->radius = best_radius;
    st_node->error_best_cut = min_error;
}

void SliceTreeSamp::optimal_cut_exact(st_node_t* st_node) const
{
    unsigned int best_center = st_node->partition[0];
    unsigned int best_radius = 0;
    double min_error = std::numeric_limits<double>::max();
    std::pair<double , unsigned int> e_r;
    double error;
    unsigned int radius;

    /*Computes the best radius for each possible center in the
    * partition*/
    for(unsigned int c = 0; c < st_node->partition.size(); c++)
    {
        e_r = min_error_radius(st_node->partition[c],
            st_node->partition, st_node->diameter,
            st_node->in_partition, st_node->average);
        error = global_error - e_r.first;
        radius = e_r.second;

        if(error < min_error)
        {
            best_center = st_node->partition[c];
            best_radius = radius;
            min_error = error;
        }
    }

    st_node->center = best_center;
    st_node->radius = best_radius;
    st_node->error_best_cut = min_error;
}

/**
 * Computes theta, which is the range in which all values are.
 * @param
 * @return
 * @throws
**/
double SliceTreeSamp::compute_theta()
{
    double min_value = std::numeric_limits<double>::max();
    double max_value = -1*std::numeric_limits<double>::max();

    for(unsigned int v = 0; v < graph->size(); v++)
    {
        if(graph->orig_value(v) < min_value)
        {
            min_value = graph->orig_value(v);
        }

        if(graph->orig_value(v) > max_value)
        {
            max_value = graph->orig_value(v);
        }
    }

    theta = fabs(max_value - min_value);

    return theta;
}

/**
 * Computes a lower bound on the size of a partition. Computing the actual
 * number of vertices inside and outside a slice can be expensive and there
 * is no way we can keep this information for slices other than the first
 * one. Therefore, we use this simple bound that returns the size of the partition
 * for the largest first slice considering a radius that is smaller or equal to the
 * radius of the slice of interest but that cannot intersect with any existing slice.
 * @param center center
 * @param radius radius
 * @param partition partition
 * @return lower bound on the size of the partition defined by the given
 * center and radius.
 * @throws
**/
const unsigned int SliceTreeSamp::lower_bound_size_partition(
    const unsigned int center,
    const unsigned int radius,
    const std::vector<unsigned int>& partition) const
{
    unsigned int r = radius;

    if(dist_near_center.at(center) < std::numeric_limits<unsigned int>::max()
        && dist_near_center.at(center) - radius_near_center.at(center) < radius)
    {
        r = dist_near_center.at(center) - radius_near_center.at(center);
    }

    if(dist_center_part.at(center) < std::numeric_limits<unsigned int>::max()
        && radius_near_center.at(center) - dist_center_part.at(center) < r)
    {
        r = radius_near_center.at(center) - dist_center_part.at(center);
    }

    return graph->get_partition_size(center, r);
}

/**
 * Computes an upper bound on the size of the partition, which, basically,
 * does not consider any intersection between slices. The value is exact only
 * for the first slice.
 * @param center center
 * @param radius radius
 * @param partition partition
 * @throws
 * @return
**/
const unsigned int SliceTreeSamp::upper_bound_size_partition(
    const unsigned int center,
    const unsigned int radius,
    const std::vector<unsigned int>& partition) const
{
    unsigned int val_graph = graph->get_partition_size(center, radius);

    if(val_graph < partition.size())
    {
        return val_graph;
    }
    else
    {
        return partition.size();
    }
}

/**
 * Computes a lower bound on the size of the complement of
 * a partition, which, basically,
 * does not consider any intersection between slices. The value is exact only
 * for the first slice.
 * @param center center
 * @param radius radius
 * @param partition partition
 * @throws
 * @return
**/
const unsigned int SliceTreeSamp::lower_bound_size_comp_partition(
    const unsigned int center,
    const unsigned int radius,
    const std::vector<unsigned int>& partition) const
{
    unsigned int val_graph = graph->get_partition_size(center, radius);

    if(val_graph > partition.size())
    {
        return 1;
    }
    else
    {
        return partition.size() - val_graph;
    }
}

/**
 * Computes a lower bound on the size of the complement of
 * a partition, which, basically,
 * does not consider any intersection between slices. The value is exact only
 * for the first slice.
 * @param center center
 * @param radius radius
 * @param partition partition
 * @throws
 * @return
**/
const unsigned int SliceTreeSamp::upper_bound_size_comp_partition(
    const unsigned int center,
    const unsigned int radius,
    const std::vector<unsigned int>& partition) const
{
    unsigned int val_graph = lower_bound_size_partition(center, radius, partition);

    if(val_graph > partition.size())
    {
        return 1;
    }
    else
    {
        return partition.size() - val_graph;
    }
}

/**
 * Computes a probabilistic upper bound on the error reduction of
 * a slice based on an estimate for the mean value inside the partition,
 * which is computed using the sample.
 * @param center center
 * @param radius radius
 * @param partition partition
 * @param average partition average value
 * @param weighted_mean weighted mean of the partition generated by
 * the slice computed using the sample.
 * @param num_samples_part number of samples used to compute the
 * weighted mean.
 * @throws
 * @return upper bound.
**/
std::pair<double, double>
SliceTreeSamp::upper_bound_error_reduction_mean_estimate_in(
    const unsigned int center, const unsigned int radius,
    const std::vector<unsigned int>& partition,
    const double average, const double weighted_mean,
    const unsigned int num_samples_part) const
{
    std::pair<double,double> res;
    double epsilon = compute_epsilon_mean(theta, num_samples_part, delta);
    unsigned int size_partition = upper_bound_size_partition
        (center, radius, partition);
    unsigned int size_comp_partition = lower_bound_size_comp_partition
        (center, radius, partition);

    double bound_one = (double) (pow(average - weighted_mean
        + epsilon, 2) * size_partition * partition.size())
        / size_comp_partition;
    double bound_two = (double) (pow(average - weighted_mean
        - epsilon, 2) * size_partition * partition.size())
        / size_comp_partition;

    double estimate = 0;

    if(bound_one > bound_two)
    {
//      printf("**center = %d, radius = %d, bound = %lf, weighted_mean = %lf, epsilon = %lf, num_samples = %d, average = %lf, size_partition = %d, size_comp_partition = %d, partition.size = %d\n",
//          center, radius, bound_one, weighted_mean, epsilon, num_samples_part, average, size_partition, size_comp_partition, partition.size());
        res.first = bound_one;
    }
    else
    {
//      printf("**center = %d, radius = %d, bound = %lf, weighted_mean = %lf, epsilon = %lf, num_samples = %d, average = %lf, size_partition = %d, size_comp_partition = %d, partition.size = %d\n",
//          center, radius, bound_two, weighted_mean, epsilon, num_samples_part, average, size_partition, size_comp_partition, partition.size());
        res.first = bound_two;
    }


    if(size_comp_partition > 0)
    {
        estimate = (double) (pow(average - weighted_mean, 2)
            * size_partition * partition.size())
            / size_comp_partition;
    }

    res.second = estimate;

    return res;
}

/**
 * Computes a probabilistic upper bound on the error reduction of
 * a slice based on an estimate for the mean value outside the partition,
 * which is computed using the sample.
 * @param center center
 * @param radius radius
 * @param partition partition
 * @param average partition average value
 * @param weighted_mean weighted mean of the partition generated by
 * the slice computed using the sample.
 * @param num_samples_part number of samples used to compute the
 * weighted mean.
 * @throws
 * @return upper bound.
**/

std::pair<double,double>
    SliceTreeSamp::upper_bound_error_reduction_mean_estimate_out(
    const unsigned int center, const unsigned int radius,
    const std::vector<unsigned int>& partition,
    const double average, const double weighted_mean,
    const unsigned int num_samples_part) const
{
    std::pair<double,double> res;
    double epsilon = compute_epsilon_mean(theta, num_samples_part, delta);
    unsigned int size_partition = lower_bound_size_partition
        (center, radius, partition);
    unsigned int size_comp_partition = upper_bound_size_comp_partition
        (center, radius, partition);

    double bound_one = (double) (pow(average - weighted_mean
        + epsilon, 2) * size_comp_partition * partition.size())
        / size_partition;
    double bound_two = (double) (pow(average - weighted_mean
        - epsilon, 2) * size_comp_partition * partition.size())
        / size_partition;


    if(bound_one > bound_two)
    {
//      printf("@@center = %d, radius = %d, bound = %lf, weighted_mean = %lf, epsilon = %lf, num_samples = %d, average = %lf, size_partition = %d, size_comp_partition = %d, partition.size = %d\n",
//          center, radius, bound_one, weighted_mean, epsilon, num_samples_part, average, size_partition, size_comp_partition, partition.size());
        res.first = bound_one;
    }
    else
    {
//      printf("@@center = %d, radius = %d, bound = %lf, weighted_mean = %lf, epsilon = %lf, num_samples = %d, average = %lf, size_partition = %d, size_comp_partition = %d, partition.size = %d\n",
//          center, radius, bound_two, weighted_mean, epsilon, num_samples_part, average, size_partition, size_comp_partition, partition.size());
        res.first = bound_two;
    }

    double estimate = (double) (pow(average - weighted_mean, 2)
        * size_comp_partition * partition.size())
        / size_partition;
    res.second = estimate;

    return res;
}

/**
 * Computes a probabilistic upper bound on the error reduction of
 * a slice based on the number of vertices inside the partition sampled
 * in a biased sample.
 * @param center center
 * @param radius radius
 * @param partition partition
 * @param num_samples_part number of samples used to compute the
 * weighted mean.
 * @throws
 * @return upper bound.
**/
std::pair<double,double>
    SliceTreeBiasSamp::upper_bound_error_reduction_num_samples
    (const unsigned int center, const unsigned int radius,
    const std::vector<unsigned int>& partition,
    const unsigned int num_samples_part,
    const unsigned int total_samples) const
{
    std::pair<double,double> res;
    double estimate = 0;
    double sampling_rate = (double) num_samples_part / total_samples;
    double epsilon = compute_epsilon_num_samples(total_samples, delta);
    unsigned int size_partition = lower_bound_size_partition
        (center, radius, partition);
    unsigned int size_comp_partition = lower_bound_size_comp_partition
        (center, radius, partition);

    double bound = (double) (pow(sampling_rate + epsilon, 2)
        * pow(graph->get_lambda(), 2) * partition.size())
        / (size_partition * size_comp_partition);

    if(size_comp_partition > 0)
    {
        estimate = (double) (pow(sampling_rate, 2)
            * pow(graph->get_lambda(), 2) * partition.size())
            / (size_partition * size_comp_partition);
    }

    res.first = bound;
    res.second = estimate;
//  printf("##center = %d, radius = %d, bound = %lf, rate = %lf, epsilon = %lf, size_partition = %d, size_complement_partition = %d, lambda = %lf, num_samples_part = %d, num_samples = %d, size_x = %d\n",
//      center, radius, bound, sampling_rate, epsilon, size_partition, size_comp_partition, graph->get_lambda(), num_samples_part, num_samples, partition.size());

    return res;
}

/* Computes a probabilistic upper bound on the error reduction of
 * a slice based on the number of vertices inside the partition sampled
 * in a biased sample. Because this is not biased sampling, just returns
 * a very large number.
 * @param center center
 * @param radius radius
 * @param partition partition
 * @param num_samples_part number of samples used to compute the
 * weighted mean.
 * @throws
 * @return upper bound.
**/
std::pair<double, double>
    SliceTreeUnifSamp::upper_bound_error_reduction_num_samples
    (const unsigned int center, const unsigned int radius,
    const std::vector<unsigned int>& partition,
    const unsigned int num_samples_part,
    const unsigned int total_samples) const
{
    std::pair<double, double> res;
    res.first = std::numeric_limits<double>::max();
    res.second = 0;
    return res;
}

double SliceTreeBiasSamp::compute_estimate(const double one,
    const double two, const double three) const
{
    double estimate = (double) 3 /
        (((double) 1 / one)
        +((double) 1 / two)+
        (double) 1 / three);

    return estimate;
}

double SliceTreeUnifSamp::compute_estimate(const double one,
    const double two, const double three) const
{
    double estimate = (double) 2 /
        (((double) 1 / one)
        +((double) 1 / two));

    return estimate;
}

/**
 * Computes upper bounds on the error reduction of slices centered
 * at a given vertex using sampling and inserts them into a set of
 * upper bounds.
 * @param upper_bounds set of upper bounds
 * @param center center
 * @param partition partition
 * @param diameter diameter
 * @param in_partition bitmap for the partition
 * @param average partition average
 * @return
 * @throws
**/
void SliceTreeSamp::upper_bound_error_reduction(up_bound_t* up_bound,
    const unsigned int center,
    const std::vector<unsigned int>& partition,
    const unsigned int diameter,
    const std::vector<bool>& in_partition,
    const double average,
    const double sum_weighted_values,
    const double sum_weights,
    const unsigned int total_samples,
    const double sse_partition) const
{
    std::list<unsigned int>* vertices_at_dist_r;
    double sum_weights_in = 0;
    double sum_weighted_values_in = 0;
    double sum_weights_out;
    double sum_weighted_values_out;
    unsigned int num_samples_part_in = 0;
    unsigned int num_samples_part_out;
    double weighted_mean_in = 0;
    double weighted_mean_out = 0;
    std::pair<double, double> up_est_one;
    std::pair<double, double> up_est_two;
    std::pair<double, double> up_est_three;
    unsigned int vertex;
    double bound_value;
    unsigned int max_radius_slice;
    double estimate;

    if(max_radius > diameter)
    {
        max_radius_slice = diameter;
    }
    else
    {
        max_radius_slice = max_radius;
    }

    /*For radius = 0, we compute the actual reduction, instead of an
     * upper bound.*/
    up_bound->center = center;

    if(partition.size() - 1 > 0)
    {
        up_bound->bound =
            (double)(pow(graph->orig_value(partition.at(center)) - average, 2) * partition.size()) / (partition.size() - 1);
        up_bound->estimate = up_bound->bound;
    }
    else
    {
        up_bound->bound = 0;
        up_bound->estimate = 0;
    }

    up_bound->bounds.push_back(new std::pair<unsigned int, double>(0, up_bound->bound));

//  printf("center = %d, max_radius = %d, diameter = %d\n", partition.at(center), max_radius_slice, diameter);


    for(unsigned int r = 0; r <= max_radius_slice; r++)
    {
        if(r < graph->max_distance(partition.at(center)))
        {
            up_bound->radius = r;
            vertices_at_dist_r = graph->vertices_at_distance
                (partition.at(center), r);

            for(std::list<unsigned int>::iterator it =
                vertices_at_dist_r->begin();
                it != vertices_at_dist_r->end();++it)
            {
                vertex = *it;

                if(in_partition[vertex])
                {
                    sum_weighted_values_in +=
                        graph->count(vertex) * graph->value(vertex);

                    sum_weights_in +=
                        graph->count(vertex) * graph->weight(vertex);
//                  printf("+v = %d, count = %d\n", vertex, graph->count(vertex));
                    num_samples_part_in += graph->count(vertex);
                }
            }
        }

        if(sum_weights_in > 0)
        {
            weighted_mean_in = (double) sum_weighted_values_in
                / sum_weights_in;
        }

        sum_weights_out = sum_weights - sum_weights_in;
        sum_weighted_values_out =
            sum_weighted_values - sum_weighted_values_in;

        if(sum_weights_out > 0)
        {
            weighted_mean_out = (double) sum_weighted_values_out
                / sum_weights_out;
        }
        else
        {
            weighted_mean_out = 0;
        }

        num_samples_part_out = total_samples - num_samples_part_in;

        up_est_one =
            upper_bound_error_reduction_mean_estimate_in
                (partition.at(center), r,
                partition, average,
                weighted_mean_in, num_samples_part_in);

        up_est_two =
            upper_bound_error_reduction_mean_estimate_out
                (partition.at(center), r,
                partition, average,
                weighted_mean_out, num_samples_part_out);

        up_est_three = upper_bound_error_reduction_num_samples
            (partition.at(center), r, partition, num_samples_part_in,
            total_samples);

        if(r > 0)
        {
            if(up_est_one.first < up_est_two.first)
            {
                if(up_est_one.first < up_est_three.first)
                {
                    if(up_est_one.first < sse_partition)
                    {
                        bound_value = up_est_one.first;
                        num_pruned_bound_1++;
                    }
                    else
                    {
                        bound_value = sse_partition;
                    }
                }
                else
                {
                    if(up_est_three.first < sse_partition)
                    {
                        bound_value = up_est_three.first;

                        if(up_est_three.first < up_est_one.first)
                        {
                            num_pruned_bound_3++;
                        }
                    }
                    else
                    {
                        bound_value = sse_partition;
                    }
                }
            }
            else
            {
                if(up_est_two.first < up_est_three.first)
                {
                    if(up_est_two.first < sse_partition)
                    {
                        bound_value = up_est_two.first;
                        num_pruned_bound_2++;
                    }
                    else
                    {
                        bound_value = sse_partition;
                    }
                }
                else
                {
                    if(up_est_three.first < sse_partition)
                    {
                        bound_value = up_est_three.first;

                        if(up_est_three.first < up_est_two.first)
                        {
                            num_pruned_bound_3++;
                        }
                    }
                    else
                    {
                        bound_value = sse_partition;
                    }
                }
            }

            up_bound->bounds.push_back(
                new std::pair<unsigned int, double>
                    (r, bound_value));

            if(bound_value > up_bound->bound)
            {
                up_bound->bound = bound_value;
            }

            estimate = compute_estimate(up_est_one.second,
                up_est_two.second, up_est_three.second);

//          printf("estimate=%lf, one=%lf, second=%lf, third=%lf\n", estimate,
//              up_est_one.second, up_est_two.second, up_est_three.second);
//          printf("bound=%lf, one=%lf, second=%lf, third=%lf\n", bound_value,
//              up_est_one.first, up_est_two.first, up_est_three.first);

            if(estimate > up_bound->estimate)
            {
                up_bound->estimate = estimate;
            }
        }
    }

//  printf("--radius = %d bound = %lf\n", up_bound->radius, up_bound->bound);
}

void free_upper_bound(up_bound_t* up_bound)
{
    std::list< std::pair<unsigned int, double> *>::iterator it;
    for(it = up_bound->bounds.begin();
        it != up_bound->bounds.end(); ++it)
    {
        delete *it;
    }

    delete up_bound;
}

void SliceTreeSamp::optimal_cut(st_node_t* st_node)
{
    std::list<up_bound_t*> upper_bounds;
    up_bound_t* up_bound;
    up_bound_t* up_bound_best_estimate;
    std::list<up_bound_t*>::iterator u;
    double best_estimate = -1*std::numeric_limits<double>::max();
    std::pair<double , unsigned int> e_r;
    double sum_weights;
    double sum_weighted_values;
    unsigned int total_samples = 0;
    unsigned int _total_slices = 0;
    unsigned int round = 0;
    unsigned int remaining_centers = st_node->partition.size();
    double opt_reduction = -1*std::numeric_limits<double>::max();
    unsigned int opt_center = st_node->partition.at(0);
    unsigned int opt_radius = 0;
    std::list< std::pair<unsigned int, double> *>::iterator it;
    unsigned int radius;
    unsigned int center;
    unsigned int _num_pruned = 0;
    unsigned int num_samples = ceil(sampling_rate * st_node->partition.size());

    while(remaining_centers > num_samples)
    {
        //printf("remaining centers = %d, samples = %d\n", remaining_centers, total_samples);
        best_estimate = -1*std::numeric_limits<double>::max();
        up_bound_best_estimate = NULL;

        if(round == 0)
        {
            graph->set_sample(num_samples, st_node->partition);
            total_samples += num_samples;
            theta = graph->get_theta();
            graph->build_distance_str_slice_tree_sample(
                max_radius, st_node->partition);
            sum_weights = graph->get_sum_weights();
            sum_weighted_values = graph->get_sum_weighted_values();

            /*Computing upper and lower bounds on the error reduction for
             * all slices in partition using sampling*/
            for(unsigned int c = 0; c < st_node->partition.size(); c++)
            {
                up_bound = new up_bound_t;
                upper_bound_error_reduction(up_bound, c,
                    st_node->partition, st_node->diameter,
                    st_node->in_partition, st_node->average,
                    sum_weighted_values, sum_weights,
                    total_samples, st_node->error_partition);
                upper_bounds.push_back(up_bound);
                _total_slices += upper_bounds.back()->bounds.size();
//              printf("center:%d,estimate = %lf\n", c, up_bound->estimate);
//              printf("center:%d,upper_bound:%lf,radius:%d\n", c, up_bound->bound, up_bound->radius);
                if(up_bound->estimate >= best_estimate)
                {
                    best_estimate = up_bound->estimate;
                    up_bound_best_estimate = up_bound;
                }
            }

            _num_pruned = _total_slices;
        }
        else
        {
            graph->resample(num_samples, st_node->partition);
            total_samples += num_samples;
            graph->build_distance_str_slice_tree_sample(
                max_radius, st_node->partition);
            sum_weights = graph->get_sum_weights();
            sum_weighted_values = graph->get_sum_weighted_values();

            for(u = upper_bounds.begin(); u != upper_bounds.end(); u++)
            {
                up_bound = *u;
                center = up_bound->center;
                radius = up_bound->radius;
                free_upper_bound(up_bound);
                up_bound = new up_bound_t;
                *u = up_bound;
                upper_bound_error_reduction(up_bound, center,
                    st_node->partition, radius,
                    st_node->in_partition, st_node->average,
                    sum_weighted_values, sum_weights,
                    total_samples, st_node->error_partition);

//              printf("center:%d,estimate = %lf\n", center, up_bound->estimate);
//              printf("center:%d,upper_bound:%lf,radius:%d\n", center, up_bound->bound, up_bound->radius);
                if(up_bound->estimate >= best_estimate)
                {
                    up_bound_best_estimate = up_bound;
                    best_estimate = up_bound->estimate;
                }
            }
        }

        //printf("opt_reduction = %lf, best_estimate = %lf\n", opt_reduction, best_estimate);

        if(best_estimate >= opt_reduction)
        {
            _num_pruned -= up_bound_best_estimate->radius;

            /*Computing the actual error reduction*/
            e_r = min_error_radius(
                st_node->partition.at(up_bound_best_estimate->center),
                st_node->partition, up_bound_best_estimate->radius,
                st_node->in_partition, st_node->average);

//          printf("computed: c=%d,r=%d,est=%lf,red=%lf\n",
//              st_node->partition.at(up_bound_best_estimate->center),
//              up_bound_best_estimate->radius, best_estimate, e_r.first);

            if(e_r.first > opt_reduction)
            {
                opt_reduction = e_r.first;
                opt_radius = e_r.second;
                opt_center =
                    st_node->partition.at(
                        up_bound_best_estimate->center);
            }

            up_bound_best_estimate->bound = 0;
        }

        /*Pruning*/
        u = upper_bounds.begin();

        while(u != upper_bounds.end())
        {
            up_bound = *u;

            if(rho * up_bound->bound <= opt_reduction)
            {
                free_upper_bound(up_bound);
                u = upper_bounds.erase(u);
                remaining_centers--;
            }
            else
            {
                radius = 0;

                for(it = up_bound->bounds.begin();
                    it != up_bound->bounds.end();++it)
                {
                    if((*it)->second * rho > opt_reduction)
                    {
                        if(radius < (*it)->first)
                        {
                            radius = (*it)->first;
                        }
                    }
                }

                up_bound->radius = radius;

                ++u;
            }
        }

        round++;
    }

//  printf("Computing actual reductions...\n");

    u = upper_bounds.begin();

    while(u != upper_bounds.end())
    {
        up_bound = *u;

        if(rho * up_bound->bound > opt_reduction)
        {
            radius = 0;

            for(it = up_bound->bounds.begin();
                it != up_bound->bounds.end();++it)
            {
                if((*it)->second * rho > opt_reduction)
                {
                    if(radius < (*it)->first)
                    {
                        radius = (*it)->first;
                    }
                }
            }

            _num_pruned -= radius;

            e_r = min_error_radius(
                st_node->partition.at(up_bound->center),
                st_node->partition, radius,
                st_node->in_partition, st_node->average);

//          printf("computed: c=%d,r=%d,red=%lf\n",
//              st_node->partition.at(up_bound->center),
//              up_bound->radius, e_r.first);


            if(e_r.first > opt_reduction)
            {
                opt_reduction = e_r.first;
                opt_radius = e_r.second;
                opt_center =
                    st_node->partition.at(
                        up_bound->center);
            }
        }

        free_upper_bound(up_bound);
        u = upper_bounds.erase(u);
    }

    st_node->center = opt_center;
    st_node->radius = opt_radius;
    st_node->error_best_cut = global_error - opt_reduction;
    num_pruned += _num_pruned;
    total_slices += _total_slices;
    //printf("pruning rate = %lf\n", double(num_pruned)/total_slices);

//  printf("**center = %d, radius = %d, reduction = %lf\n", opt_center, opt_radius, opt_reduction);
    //printf("$$best: c=%d,r=%d,red=%lf\n", opt_center, opt_radius, opt_reduction);
}

/**
 * Computes the difference coefficients for a wavelet-like
 * decomposition of the slice tree for a slice tree node
 * @param st_node slice tree node
 * @return
 * @throws
**/
void compute_difference_coefficients_st_node(st_node_t* st_node)
{
    if(st_node->left != NULL and st_node->right != NULL)
    {
        st_node->difference = (st_node->average - st_node->left->average) * st_node->left->size;

        compute_difference_coefficients_st_node(st_node->left);
        compute_difference_coefficients_st_node(st_node->right);
    }
}

/**
 * Computes the difference coefficients for a wavelet-like
 * decomposition of the slice tree
 * @param
 * @return
 * @throws
**/
void SliceTree::compute_difference_coefficients()
{
    compute_difference_coefficients_st_node(tree);
}

/**
 * Identifies the optimal radius for a given center and partition
 * @param center center to be considered
 * @param partition to be split
 * @param diameter of the partition to be split
 * @param in_partition bitmap of the partition to be split
 * @return pair <error, radius>
 * @throws
**/
const std::pair<double, unsigned int> SliceTree::min_error_radius(const unsigned int center,
    const std::vector<unsigned int>& partition, unsigned int max_radius_slice,
    const std::vector<bool>& in_partition, const double average) const
{
    double curr_sum_values = 0;
    double curr_average = 0;
    unsigned int vertex;
    std::list<unsigned int>* vertices_at_dist_r;

    double max_reduction = 0;
    double reduction;
    unsigned int radius = 0;

    unsigned int curr_num_vertices = 0;
    std::vector<std::list<unsigned int>*> vertices_at_distance;

    if(max_radius_slice < max_radius)
    {
        graph->build_distance_str_slice_tree_vertex(center, vertices_at_distance, max_radius_slice);
    }
    else
    {
        graph->build_distance_str_slice_tree_vertex(center, vertices_at_distance, max_radius);
    }

    /*Computing the sse for each possible radius*/
    for(unsigned int r = 0; r < vertices_at_distance.size(); r++)
    {
        /*Getting the list of vertices at distance r from the center*/
        vertices_at_dist_r = vertices_at_distance[r];

        for(std::list<unsigned int>::iterator it = vertices_at_dist_r->begin();
            it != vertices_at_dist_r->end();++it)
        {
            vertex = *it;

            if(in_partition[vertex])
            {
                curr_sum_values = curr_sum_values
                    + graph->orig_value(vertex);
                curr_num_vertices = curr_num_vertices + 1;
            }
        }

        if(curr_num_vertices > 0)
        {
            curr_average = (double) curr_sum_values / curr_num_vertices;
        }
        else
        {
            curr_average = 0;
        }

        if(curr_num_vertices < partition.size())
        {
            reduction = (double) (pow(curr_average - average, 2)
                * curr_num_vertices * partition.size())
                / (partition.size() - curr_num_vertices);
        }
        else
        {
            reduction = 0;
        }

        if(reduction > max_reduction)
        {
            max_reduction = reduction;
            radius = r;
        }
    }

    for(unsigned int v = 0; v < vertices_at_distance.size(); v++)
    {
        delete vertices_at_distance.at(v);
    }

    std::pair<double, unsigned int> error_radius;
    error_radius.first = max_reduction;
    error_radius.second = radius;

    return error_radius;
}

/**
 * Wrapper for compress heuristics
 * @param
 * @return
 * @throws
**/
void SliceTree::compress(const unsigned int budget){
    if (exhaustive_split) compressExhaustive(budget);
    else compressGreedy(budget);
}

/**
 * Runs the slice tree compression. All partitions probed
 * @param
 * @return
 * @throws
**/
void SliceTree::compressExhaustive(const unsigned int budget)
{
    budget_compression = budget;
    tree = new st_node_t;
    n_partitions = num_partitions(budget_compression,
        graph->size(), graph->diameter());

    /*Setting the initial partition, which contains all the vertices
    * in the graph*/
    tree->partition.reserve(graph->size());
    tree->in_partition.reserve(graph->size());

    for(unsigned int v = 0; v < graph->size(); v++)
    {
        tree->partition.push_back(v);
        tree->in_partition.push_back(true);
    }

    tree->size = tree->partition.size();
    tree->center = 0;
    tree->radius = UINT_MAX;
    tree->difference = 0;
    tree->left = NULL;
    tree->right = NULL;
    tree->error_partition = sse_partition(tree->partition);
    tree->average = average_partition(tree->partition);
    tree->diameter = graph->diameter();

    /*The initial global error is the error of the initial partition*/
    global_error = tree->error_partition;

    /*Computing the best cut for the initial partition*/
    optimal_cut(tree);

    /*Candidate cuts is a min-heap, we use it to find the next best cut*/
    std::vector<st_node_t*> candidate_cuts;
    candidate_cuts.reserve(n_partitions);
    candidate_cuts.push_back(tree);

    std::push_heap(candidate_cuts.begin(), candidate_cuts.end(), CompareCuts());

    unsigned int n_part = 1;
    st_node_t* part;

    while(n_part < n_partitions && candidate_cuts.size() > 0)
    {
        part = candidate_cuts.front();
        std::pop_heap (candidate_cuts.begin(), candidate_cuts.end());
        candidate_cuts.pop_back();

        /*Performing the best cut, in case the cut is trivial
        (i.e. it generates an empty partition) the split is not performed*/
        if(split_partition(part))
        {
            /*Processing new partitions*/

            if(n_part + 1 < n_partitions)
            {
                optimal_cut(part->left);
                candidate_cuts.push_back(part->left);

                optimal_cut(part->right);
                candidate_cuts.push_back(part->right);

                std::push_heap(candidate_cuts.begin(), candidate_cuts.end(), CompareCuts());
            }

            n_part = n_part + 1;
        }
    }

    compute_difference_coefficients();
}


/**
 * Runs the slice tree compression. Only the highest error slice is
 * probed for all possible cuts
 * @param
 * @return
 * @throws
**/
void SliceTree::compressGreedy(const unsigned int budget)
{
    budget_compression = budget;
    tree = new st_node_t;
    n_partitions = num_partitions(budget_compression,
        graph->size(), graph->diameter());

    /*Setting the initial partition, which contains all the vertices
    * in the graph*/
    tree->partition.reserve(graph->size());
    tree->in_partition.reserve(graph->size());

    for(unsigned int v = 0; v < graph->size(); v++)
    {
        tree->partition.push_back(v);
        tree->in_partition.push_back(true);
    }

    tree->size = tree->partition.size();
    tree->center = 0;
    tree->radius = UINT_MAX;
    tree->difference = 0;
    tree->left = NULL;
    tree->right = NULL;
    tree->error_partition = sse_partition(tree->partition);
    tree->average = average_partition(tree->partition);
    tree->diameter = graph->diameter();

    /*The initial global error is the error of the initial partition*/
    global_error = tree->error_partition;

    /*Candidate cuts is a min-heap, we use it to find the next best cut*/
    std::vector<st_node_t*> candidate_cuts;
    candidate_cuts.reserve(n_partitions);


    unsigned int n_part = 1;
    st_node_t* part;
    st_node_t* toCut = getMaxSSEPartiton(tree);


    while(n_part < n_partitions && toCut!=NULL)
    {
        /*Computing the best cut for the initial partition*/
        optimal_cut(toCut);
        candidate_cuts.push_back(toCut);
        std::push_heap(candidate_cuts.begin(), candidate_cuts.end(), CompareCuts());

        part = candidate_cuts.front();
        std::pop_heap (candidate_cuts.begin(), candidate_cuts.end());
        candidate_cuts.pop_back();

        /*Performing the best cut, in case the cut is trivial
        (i.e. it generates an empty partition) the split is not performed*/
        if(split_partition(part))
        {
            /*Processing new partitions*/

            if(n_part + 1 < n_partitions)
            {
                // get max sse partition
                toCut = getMaxSSEPartiton(tree);
            }

            n_part = n_part + 1;
        }
    }

    compute_difference_coefficients();
}

/**
 * Find partition of max SSE
 * @param
 * @return
 * @throws
**/
st_node_t* SliceTree::getMaxSSEPartiton(st_node_t* root){
    if (root->left==NULL && root->right==NULL) {
        return root;
    } else {
        st_node_t* maxErrHere = getMaxSSEPartiton(root->left);
        st_node_t* maxErrRight = getMaxSSEPartiton(root->right);
        if (maxErrRight->error_partition > maxErrHere->error_partition)
            maxErrHere = maxErrRight;
        return maxErrHere;
    }
}

/**
 * Frees the memory of the tree recursively
 * @param st_node parent node
 * @return
 * @throws
**/
void free_st_node(st_node_t* st_node)
{
    if(st_node->left != NULL)
    {
        free_st_node(st_node->left);
    }

    if(st_node->right != NULL)
    {
        free_st_node(st_node->right);
    }

    delete st_node;
}

/**
 * Destructor
 * @param
 * @return
 * @throws
**/
SliceTree::~SliceTree()
{
    free_st_node(tree);
}

/**
 * Extends a slice tree node. A non-leaf node from slice tree recovered from a
 * serialized file is extended with left and right children.
 * @param
 * @return
 * @throws
**/
void SliceTree::extend_st_node(st_node_t* st_node, Graph* graph)
{
    bool extend_left = true;
    bool extend_right = true;

    if(st_node->left == NULL)
    {
        st_node->left = new st_node_t;
        extend_left = false;
        st_node->left->difference = 0;
        st_node->left->center = 0;
        st_node->left->radius = 0;
        st_node->left->error_partition = 0;
        st_node->left->error_best_cut = 0;
    }

    if(st_node->right == NULL)
    {
        st_node->right = new st_node_t;
        extend_right = false;
        st_node->right->difference = 0;
        st_node->right->center = 0;
        st_node->right->radius = 0;
        st_node->right->error_partition = 0;
        st_node->right->error_best_cut = 0;
    }

    /*Setting the partitions for children nodes*/
    st_node->left->in_partition.reserve(graph->size());
    st_node->right->in_partition.reserve(graph->size());
    st_node->left->partition.reserve(st_node->partition.size());
    st_node->right->partition.reserve(st_node->partition.size());

    for(unsigned int v = 0; v < graph->size(); v++)
    {
        st_node->left->in_partition.push_back(false);
        st_node->left->in_partition.push_back(false);
    }

    graph->bounded_bfs(st_node->left->partition, st_node->center, st_node->radius, st_node->in_partition);

    unsigned int vertex;

    for(unsigned int v = 0; v < st_node->left->partition.size(); v++)
    {
        vertex = st_node->left->partition[v];
        st_node->left->in_partition[vertex] = true;
    }

    for(unsigned int v = 0; v < st_node->partition.size(); v++)
    {
        vertex = st_node->partition[v];

        if(! st_node->left->in_partition[vertex])
        {
            st_node->right->in_partition[vertex] = true;
            st_node->right->partition.push_back(vertex);
        }
    }


    st_node->left->average = st_node->average
        - (double) st_node->difference / st_node->left->partition.size();

    st_node->right->average = st_node->average
        + (double) st_node->difference / st_node->right->partition.size();

    st_node->left->error_partition = sse_partition(st_node->left->partition);
    st_node->right->error_partition = sse_partition(st_node->right->partition);

    if(extend_left)
    {
        extend_st_node(st_node->left, graph);
    }
    else
    {
        st_node->left->left = NULL;
        st_node->left->right = NULL;
    }

    if(extend_right)
    {
        extend_st_node(st_node->right, graph);
    }
    else
    {
        st_node->right->left = NULL;
        st_node->right->right = NULL;
    }
}

/**
 * Extends the slice tree recovered from a serialized file
 * @param
 * @return
 * @throws
**/
void SliceTree::extend_tree()
{
    tree->partition.reserve(graph->size());
    tree->in_partition.reserve(graph->size());

    for(unsigned int v = 0; v < graph->size(); v++)
    {
        tree->partition.push_back(v);
        tree->in_partition.push_back(true);
    }

    extend_st_node(tree, graph);
}

/**
 * Clears all the partition for a st_node and does it recursively for the
 * children nodes.
 * @param
 * @return
 * @throws
**/
void clear_partitions_st_node(st_node_t* st_node)
{
    st_node->partition.clear();

    if(st_node->left != NULL)
    {
        if(st_node->left->left == NULL && st_node->left->right == NULL)
        {
            delete st_node->left;
            st_node->left = NULL;
        }
        else
        {
            clear_partitions_st_node(st_node->left);
        }
    }

    if(st_node->right != NULL)
    {
        if(st_node->right->left == NULL && st_node->right->right == NULL)
        {
            delete st_node->right;
            st_node->right = NULL;
        }
        else
        {
            clear_partitions_st_node(st_node->right);
        }
    }
}

/**
 * Clears all the partition in the slice tree
 * recursivelly.
 * @param
 * @return
 * @throws
**/
void SliceTree::clear_partitions()
{
    clear_partitions_st_node(tree);
}

/**
 * Splits a partition given the center and radius defined
 * by the slice tree node
 * @param st_node slice tree node
 * @return true in case the split was performed
 *  false, otherwise
 * @throws
**/
const bool SliceTree::split_partition(st_node_t* st_node)
{
    st_node_t* left = new st_node_t;
    st_node_t* right = new st_node_t;

    left->in_partition.reserve(graph->size());
    right->in_partition.reserve(graph->size());
    left->partition.reserve(st_node->partition.size());
    right->partition.reserve(st_node->partition.size());

    left->left = NULL;
    left->right = NULL;
    right->left = NULL;
    right->right = NULL;

    left->diameter = 2*st_node->radius;
    right->diameter = st_node->diameter;

    /*Performs a bounded bfs over the graph to find the vertices in the partition
    * that are at distance radius from the center*/
    graph->bounded_bfs(left->partition, st_node->center, st_node->radius, st_node->in_partition);

    /*In case one of the partitions generated by the split is empty
    * the split is not performed and the method returns false*/
    if(left->partition.size() == st_node->partition.size())
    {
        delete left;
        delete right;
        return false;
    }

    /*Computing lists and bitmaps for the new partitions*/
    for(unsigned int v = 0; v < graph->size(); v++)
    {
        left->in_partition.push_back(false);
        right->in_partition.push_back(false);
    }

    unsigned int vertex;

    for(unsigned int v = 0; v < left->partition.size(); v++)
    {
        vertex = left->partition[v];
        left->in_partition[vertex] = true;
    }

    for(unsigned int v = 0; v < st_node->partition.size(); v++)
    {
        vertex = st_node->partition[v];

        if(! left->in_partition[vertex])
        {
            right->in_partition[vertex] = true;
            right->partition.push_back(vertex);
        }
    }

    left->center = left->partition[0];
    right->center = right->partition[0];
    left->radius = 0;
    right->radius = 0;
    left->difference = 0;
    right->difference = 0;

    /*Computing averages and errors for the new partitions*/
    left->average = average_partition(left->partition);
    right->average = average_partition(right->partition);

    left->error_partition = sse_partition(left->partition);
    right->error_partition = sse_partition(right->partition);

    /*Updating partition sizes*/
    left->size = left->partition.size();
    right->size = right->partition.size();

    /*Updating the global error of the slice tree*/
    global_error = global_error - st_node->error_partition + left->error_partition + right->error_partition;

    /*Clearing some information from the parent slice tree node*/
    st_node->partition.clear();
    st_node->in_partition.clear();
    st_node->partition.reserve(0);
    st_node->in_partition.reserve(0);

    st_node->left = left;
    st_node->right = right;

    return true;
}

/**
 * Splits a partition given the center and radius defined
 * by the slice tree node when sampling is applied. The
 * difference from the standard version is that here we
 * update some data structures for computing upper and lower
 * bounds for sizes of partitions.
 * @param st_node slice tree node
 * @return true in case the split was performed
 *  false, otherwise
 * @throws
**/
const bool SliceTreeSamp::split_partition(st_node_t* st_node)
{
    st_node_t* left = new st_node_t;
    st_node_t* right = new st_node_t;

    left->in_partition.reserve(graph->size());
    right->in_partition.reserve(graph->size());
    left->partition.reserve(st_node->partition.size());
    right->partition.reserve(st_node->partition.size());

    left->left = NULL;
    left->right = NULL;
    right->left = NULL;
    right->right = NULL;

    left->diameter = 2*st_node->radius;
    right->diameter = st_node->diameter;

    /*Performs a bounded bfs over the graph to find the vertices in the partition
    * that are at distance radius from the center*/
    graph->bounded_bfs(left->partition, st_node->center, st_node->radius,
        st_node->in_partition);

    /*In case one of the partitions generated by the split is empty
    * the split is not performed and the method returns false*/
    if(left->partition.size() == st_node->partition.size())
    {
        delete left;
        delete right;
        return false;
    }

    /*Computing lists and bitmaps for the new partitions*/
    for(unsigned int v = 0; v < graph->size(); v++)
    {
        left->in_partition.push_back(false);
        right->in_partition.push_back(false);
    }

    unsigned int vertex;

    for(unsigned int v = 0; v < left->partition.size(); v++)
    {
        vertex = left->partition[v];
        left->in_partition[vertex] = true;
    }

    for(unsigned int v = 0; v < st_node->partition.size(); v++)
    {
        vertex = st_node->partition[v];

        if(! left->in_partition[vertex])
        {
            right->in_partition[vertex] = true;
            right->partition.push_back(vertex);
        }
    }

    /*Updating the nearest center information for upper-bounding the
    * sizes of partitions.*/
    graph->update_partition_size_structs(st_node->center, st_node->radius,
        dist_near_center, dist_center_part, radius_near_center,
        radius_part, st_node->in_partition);

//  graph->clean_distance_str_sample(left->partition, left->in_partition);
//  graph->clean_distance_str_sample(right->partition, right->in_partition);

    left->center = left->partition[0];
    right->center = right->partition[0];
    left->radius = 0;
    right->radius = 0;
    left->difference = 0;
    right->difference = 0;

    /*Computing averages and errors for the new partitions*/
    left->average = average_partition(left->partition);
    right->average = average_partition(right->partition);

    left->error_partition = sse_partition(left->partition);
    right->error_partition = sse_partition(right->partition);

    /*Updating partition sizes*/
    left->size = left->partition.size();
    right->size = right->partition.size();

    /*Updating the global error of the slice tree*/
    global_error = global_error - st_node->error_partition + left->error_partition + right->error_partition;

    /*Clearing some information from the parent slice tree node*/
    st_node->partition.clear();
    st_node->in_partition.clear();
    st_node->partition.reserve(0);
    st_node->in_partition.reserve(0);

    st_node->left = left;
    st_node->right = right;

    return true;
}

/**
 * Writes a slice tree node to a file
 * @param st_node slice tree node
 * @param input_file input file
 * @return
 * @throws
**/
void write_st_node(st_node_t* st_node, std::ofstream& output_file)
{
    if(st_node->left != NULL and st_node->right != NULL)
    {
        uint32_t center = (uint32_t) st_node->center;

        output_file.write(reinterpret_cast<char*>(&center), sizeof(uint32_t));

        if(st_node->radius > 255)
        {
            std::cerr << "Error: Radius is too large to be represented as an unsigned char!" << std::endl;
            exit(1);
        }

        unsigned char radius = (char) st_node->radius;

        output_file.write(reinterpret_cast<char*>(&radius), sizeof(char));

        output_file.write(reinterpret_cast<char*>(&(st_node->difference)), sizeof(float));

        unsigned char marker;

        if(st_node->right->left == NULL and st_node->right->right == NULL)
        {
            if(st_node->left->left == NULL and st_node->left->right == NULL)
            {
                marker = 0;
            }
            else
            {
                marker = 1;
            }
        }
        else
        {
            if(st_node->left->left == NULL and st_node->left->right == NULL)
            {
                marker = 2;
            }
            else
            {
                marker = 3;
            }
        }

        output_file.write(reinterpret_cast<char*>(&marker), sizeof(char));

        write_st_node(st_node->left, output_file);
        write_st_node(st_node->right, output_file);
    }
}

/**
 * Reads a slice tree node from a file
 * @param st_node slice tree node
 * @param input_file input file
 * @return
 * @throws
**/
void read_st_node(st_node_t* st_node, std::ifstream& input_file)
{
    uint32_t center;
    input_file.read(reinterpret_cast<char*>(&center), sizeof(uint32_t));
    st_node->center = (unsigned int) center;

    unsigned char radius;
    input_file.read(reinterpret_cast<char*>(&radius), sizeof(char));

    st_node->radius = (unsigned int) radius;

    input_file.read(reinterpret_cast<char*>(&st_node->difference), sizeof(float));

    unsigned char marker;

    input_file.read(reinterpret_cast<char*>(&marker), sizeof(char));

    if(marker == 0)
    {
        st_node->left = NULL;
        st_node->right = NULL;
    }
    else
    {
        if(marker == 1)
        {
            st_node->left = new st_node_t;
            read_st_node(st_node->left, input_file);
            st_node->right = NULL;
        }
        else
        {
            if(marker == 2)
            {
                st_node->right = new st_node_t;
                read_st_node(st_node->right, input_file);
                st_node->left = NULL;
            }
            else
            {
                st_node->left = new st_node_t;
                read_st_node(st_node->left, input_file);

                st_node->right = new st_node_t;
                read_st_node(st_node->right, input_file);
            }
        }
    }
}

/**
 * Writes the slice tree to a binary file.
 * Format of the binary file:
 * <average dataset><center_0><radius_0><difference_0><marker_0>...
 * <center_1><radius_1><difference_1><marker_1> ...
 * average and differences are floats (32/64 bits)
 * centers are unsigned ints (32/64 bits)
 * radii and markers are chars (4 bits)
 * @param output_file_name ouput file name
 * @return
 * @throws
**/
void SliceTree::write(const std::string& output_file_name) const
{
    std::ofstream output_file(output_file_name.c_str(), std::ios::out | std::ios::binary);

    /*Writes the overall average to the file*/
    output_file.write(reinterpret_cast<char*>(&(tree->average)), sizeof(float));

    write_st_node(tree, output_file);

    output_file.close();
}

/**
 * Buids a slice tree from the serialized content of a file
 * Format of the binary file:
 * <average dataset><center_0><radius_0><difference_0><marker_0>...
 * <center_1><radius_1><difference_1><marker_1> ...
 * average and differences are floats (32/64 bits)
 * centers are unsigned ints (32/64 bits)
 * radii and markers are chars (4 bits)
 * @param input_file_name input file
 * @return
 * @throws
**/
void SliceTree::decompress()
{
    tree = new st_node_t;
    std::ifstream input_file(compressed_file_name.c_str(), std::ios::in | std::ios::binary);

    input_file.read(reinterpret_cast<char*>(&tree->average), sizeof(float));

    read_st_node(tree, input_file);

    input_file.close();

    extend_tree();
}

/**
 * Sets values to the graph based on a node of a recovered slice tree
 * @param st_node slice tree node
 * @param graph graph
 * @return
 * @throws
**/
void set_values_graph_st_node(st_node_t* st_node, Graph* graph)
{
    unsigned int vertex;

    /*Values are set according to leaf nodes*/
    if(st_node->left == NULL && st_node->right == NULL)
    {
        for(unsigned int v = 0; v < st_node->partition.size(); v++)
        {
            vertex = st_node->partition[v];
            graph->set_value(vertex, st_node->average);
        }
    }
    else
    {
        set_values_graph_st_node(st_node->left, graph);
        set_values_graph_st_node(st_node->right, graph);
    }
}

/**
 * Computes the average value of a partition
 * @param partition partition
 * @return average
 * @throws
**/
const double SliceTree::average_partition(const std::vector<unsigned int>& partition) const
{
    double sum = 0;

    for(unsigned int v = 0; v < partition.size(); v++)
    {
        sum = sum + graph->orig_value(partition[v]);
    }

    return (double) sum / partition.size();
}

/**
 * Computes the sse of a partition
 * @param partition partition
 * @return sse
 * @throws
**/
const double SliceTree::sse_partition(const std::vector<unsigned int>& partition) const
{
    double avg = average_partition(partition);

    double sse = 0;

    for(unsigned int v = 0; v < partition.size(); v++)
    {
        sse = sse + pow(avg - graph->orig_value(partition[v]), 2);
    }


    return sse;
}

/**
 * Prints a slice tree node recursively
 * @param
 * @return
 * @throws
**/
void print_st_node(st_node_t* st_node, unsigned int depth, std::string pid,
    std::string type, Graph* graph)
{
    std::stringstream mid;
    if(st_node->left == NULL && st_node->right == NULL)
    {
        mid<<pid<<type;
    }
    else
    {
        mid<<"("<<graph->name(st_node->center)<<","
            <<st_node->radius<<","<<depth<<")" <<type;
    }

    std::cout <<pid << " " << mid.str() << " center="
        << graph->name(st_node->center) << " radius="
        << st_node->radius << " error_p="
        << st_node->error_partition << " mean=" << st_node->average;

    if(st_node->left != NULL)
    {
        std::cout<<std::endl;
        for(unsigned int d = 0; d <= depth; d++)
        {
            std::cout << " ";
        }

        print_st_node(st_node->left, depth+1, mid.str(), "i", graph);
    }

    if(st_node->right != NULL)
    {
        std::cout<<std::endl;
        for(unsigned int d = 0; d <= depth; d++)
        {
            std::cout << " ";
        }

        print_st_node(st_node->right, depth+1, mid.str(), "o", graph);
    }

    if(st_node->left == NULL && st_node->right == NULL)
    {
        std::cout<< " " ;
        for(unsigned int v = 0; v < st_node->partition.size(); v++)
        {
            std::cout << graph->name(st_node->partition[v]) << " ";
        }

        std::cout << std::endl;
    }
}

/**
 * Prints the slice tree on the terminal
 * @param
 * @return
 * @throws
**/
void SliceTree::print() const
{
    print_st_node(tree, 0, "(root)", "", graph);
}

/**
 * Sets the values (recursively) as they would be recovered after compression
 * using slice tree.
 * @param st_node slice tree node
 * @param values vector in which the values will be set
 * @return
 * @throws
**/
void set_compressed_values_st_node(st_node_t* st_node, std::vector<double>& values)
{
    unsigned int vertex;

    if(st_node->left == NULL && st_node->right == NULL)
    {
        for(unsigned int v = 0; v < st_node->partition.size(); v++)
        {
            vertex = st_node->partition[v];
            values[vertex] = st_node->average;
        }
    }
    else
    {
        set_compressed_values_st_node(st_node->left, values);
        set_compressed_values_st_node(st_node->right, values);
    }
}

/**
 * Sets the values as they would be recovered after compression
 * using slice tree.
 * @param
 * @return
 * @throws
**/
void SliceTree::set_compressed_values()
{
    values.clear();

    for(unsigned int v = 0; v < graph->size(); v++)
    {
        values.push_back(0);
    }

    clear_partitions();
    extend_tree();
    set_compressed_values_st_node(tree, values);
}

/**
 * Computes the distance between two al_nodes (partitions) as the
 * average distance between pairs of nodes from each of them.
 * @param first first average linkage node
 * @param second second average linkage node
 * @return
 * @throws
**/
double distance_al_nodes_t(al_node_t* first, al_node_t* second, Graph* graph)
{
    double average_distance = 0;
    unsigned int vertex_u;
    unsigned int vertex_v;

    for(unsigned int v = 0; v < first->partition.size(); v++)
    {
        vertex_v = first->partition.at(v);
        for(unsigned int u = 0; u < second->partition.size(); u++)
        {
            vertex_u = second->partition.at(u);
            average_distance = average_distance +
                graph->distance(vertex_v, vertex_u);
        }
    }

    average_distance = (double) average_distance /
        (first->partition.size() * second->partition.size());

    return average_distance;
}

/**
 * Computes the difference coefficients for a wavelet
 * decomposition of the average linkage tree node
 * @param al_node average linkage tree node
 * @return
 * @throws
**/
void compute_difference_coefficients_al_node_t(al_node_t* al_node)
{
    if(al_node->left != NULL and al_node->right != NULL)
    {
        al_node->difference = (al_node->average - al_node->left->average) * al_node->left->size;

        compute_difference_coefficients_al_node_t(al_node->left);
        compute_difference_coefficients_al_node_t(al_node->right);
    }
}

/**
 * Computes the average coefficients for a wavelet
 * decomposition of the average linkage tree node
 * @param al_node average linkage tree node
 * @return
 * @throws
**/
void compute_average_coefficients_al_node_t(al_node_t* al_node, Graph* graph)
{
    if(al_node->left != NULL and al_node->right != NULL)
    {
        compute_average_coefficients_al_node_t(al_node->left, graph);
        compute_average_coefficients_al_node_t(al_node->right, graph);

        al_node->average = (double)(al_node->left->size*al_node->left->average +
        al_node->right->size*al_node->right->average) /
        (al_node->left->size + al_node->right->size);
    }
    else
    {
        al_node->average = graph->value(al_node->partition.at(0));
    }
}

/**
 * Computes the difference coefficients for a wavelet
 * decomposition of the average linkage tree
 * @param
 * @return
 * @throws
**/
void AverageLinkage::compute_difference_coefficients()
{
    compute_difference_coefficients_al_node_t(tree);
}

/**
 * Computes the average coefficients for a wavelet
 * decomposition of the average linkage tree
 * @param
 * @return
 * @throws
**/
void AverageLinkage::compute_average_coefficients()
{
    compute_average_coefficients_al_node_t(tree, graph);
}

/**
 * Fills a vector with the average linkage nodes recursively
 * @param al_node average linkage node
 * @param al_nodes vector to be filled with nodes
 * @return
 * @throws
**/
void get_al_intermediate_nodes(al_node_t* al_node, std::vector<al_node_t*>& al_nodes)
{
    if(al_node->left != NULL && al_node->right != NULL)
    {
        al_nodes.push_back(al_node);
        al_node->left->average = 0;
        al_node->right->average = 0;
        get_al_intermediate_nodes(al_node->left, al_nodes);
        get_al_intermediate_nodes(al_node->right, al_nodes);
    }
}

/**
 * Wavelet coefficient pruning of the average linkage tree
 * @param
 * @return
 * @throws
**/
void AverageLinkage::keep_top_coefficients()
{
    std::vector<al_node_t*> al_nodes;
    al_nodes.reserve(graph->size());
    unsigned int num_remaining_coefficients = num_coefficients;

    /*Gets the vector of average linkage nodes*/
    get_al_intermediate_nodes(tree, al_nodes);

    /*Sorts the average linkage nodes*/
    std::sort(al_nodes.begin(), al_nodes.end(), CompareALNodes());

    /*Decides whether the global average, which is a kind of special
     * coefficient will be kept, the relevance of the average coefficient
     * is the its value times the number of leaf nodes. The relevance of
     * the other coefficients is twice their value. Don't need to multiply
     * by two if comparing only difference coefficients*/
    if(num_coefficients == 0 ||
        2*al_nodes.at(num_coefficients-1)->difference > tree->average * graph->size())
    {
        tree->average = 0;
        num_remaining_coefficients = num_coefficients;
    }
    else
    {
        num_remaining_coefficients = num_remaining_coefficients - 1;
    }

    /*Sets the coefficients not kept to 0*/
    for(unsigned int c = num_remaining_coefficients; c < al_nodes.size(); c++)
    {
        al_nodes.at(c)->difference = 0;
    }
}

/**
 * Sets the values (recursively )as they would be recovered after compression
 * using average linkage.
 * @param al_node average linkage node
 * @param values vector in which the values will be set
 * @return
 * @throws
**/
void set_compressed_values_al_node_t(al_node_t* al_node, std::vector<double>& values)
{
    if(al_node->left != NULL and al_node->right != NULL)
    {
        al_node->left->average = al_node->average
            - (double) al_node->difference / al_node->left->size;

        al_node->right->average = al_node->average
            + (double) al_node->difference / al_node->right->size;

        set_compressed_values_al_node_t(al_node->left, values);
        set_compressed_values_al_node_t(al_node->right, values);
    }
    else
    {
        values.at(al_node->partition.at(0)) = al_node->average;
    }
}

/**
 * Sets the values as they would be recovered after compression
 * using average linkage.
 * @param
 * @return
 * @throws
**/
void AverageLinkage::set_compressed_values()
{
    values.clear();

    for(unsigned int v = 0; v < graph->size(); v++)
    {
        values.push_back(0);
    }

    set_compressed_values_al_node_t(tree, values);
}

/**
 * Constructor. Builds the average linkage tree.
 * @param graph graph
 * @param budget budget
 * @return
 * @throws
**/
AverageLinkage::AverageLinkage(Graph& _graph):
    GraphCompressionAlgorithm(_graph)
{
    al_node_t* al_node;

    /*Allocating memory for the data structures*/
    partitions.reserve(graph->size());
    active_partitions.reserve(graph->size());
    distance_matrix = (double**) malloc (graph->size() * sizeof(double*));

    /*One partition is created for each vertex in the graph.
     *All partitions are set as active.*/
    for(unsigned int v = 0; v < graph->size(); v++)
    {
        al_node = new al_node_t;
        al_node->average = 0;
        al_node->difference = 0;
        al_node->size = 1;
        al_node->partition.push_back(v);
        al_node->left = NULL;
        al_node->right = NULL;

        partitions.push_back(al_node);

        distance_matrix[v] = (double*) malloc(graph->size() * sizeof(double));

        active_partitions.push_back(true);
    }

    /*Builds a distance matrix for the partitions. The diagonal is filled with infinite*/
    for(unsigned int v = 0; v < graph->size(); v++)
    {
        for(unsigned int u = 0; u < graph->size(); u++)
        {
            if(u == v)
            {
                distance_matrix[u][v] = std::numeric_limits<double>::max();
            }
            else
            {
                distance_matrix[u][v] =
                distance_al_nodes_t(partitions.at(v),
                partitions.at(u), graph);
            }
        }
    }

    /*Starts a complete descending path with the first partition (arbitrary)*/
    complete_desc_path.reserve(graph->size());
    complete_desc_path.push_back(0);

    /*Constructs a complete descending path*/
    while(construct_desc_path());

    /*The algorithm finishes when the complete descending path has size 1*/
    while(complete_desc_path.size() > 1)
    {
        /*Joins the two last partitions of the complete descending path*/
        join_partitions();

        if(! complete_desc_path.size())
        {
            for(unsigned int p = 0; p < partitions.size(); p++)
            {
                if(active_partitions.at(p))
                {
                    complete_desc_path.push_back(p);
                    break;
                }
            }
        }

        /*Constructs a complete descending path again*/
        while(construct_desc_path());
    }

    for(unsigned int v = 0; v < graph->size(); v++)
    {
        free(distance_matrix[v]);
    }

    free(distance_matrix);

    tree = partitions.at(complete_desc_path.at(0));
}

/**
 * Runs the average linkage compression.
 * @param
 * @return
 * @throws
**/
void AverageLinkage::compress(const unsigned int budget)
{
    budget_compression = budget;
    num_coefficients = (int) (floor((double) (budget_compression) /
    (floor(log2(graph->size())/8) + SIZE_FLOAT_INT)));

    if(num_coefficients > graph->size())
    {
        num_coefficients = graph->size();
    }

    compute_average_coefficients();
    compute_difference_coefficients();
    keep_top_coefficients();
}

/**
 * Tries to perform a basic operation to extend the complete descending path.
 * @param
 * @return true in case the extension should continue, false if it is stopped.
 * @throws
**/
bool AverageLinkage::construct_desc_path()
{
    /*Finds the partition that is closest to the last one in the path*/
    unsigned int last = complete_desc_path.at(complete_desc_path.size()-1);

    unsigned int closest = last;
    double distance_closest = std::numeric_limits<double>::max();

    for(unsigned int v = 0; v < graph->size(); v++)
    {
        if(active_partitions[v])
        {
            if(distance_matrix[last][v] <= distance_closest)
            {
                closest = v;
                distance_closest = distance_matrix[last][v];
            }
        }
    }

    /*Didn't find a closest partition, stop the construction*/
    if(closest == last)
    {
        return false;
    }

    /*The distance between the last and the closest is the the same as the
     *distance between the last and the partition before the last, then
     *stop the extension*/
    if(complete_desc_path.size() > 1)
    {
        unsigned int before_last =
            complete_desc_path.at(complete_desc_path.size()-2);

        if(distance_matrix[last][before_last] == distance_matrix[last][closest])
        {
            return false;
        }
    }

    /*Add the closest partition to the complete descending path*/
    complete_desc_path.push_back(closest);

    return true;
}

/**
 * Joins the two last partitions of the complete descending path.
 * @param
 * @return
 * @throws
**/
void AverageLinkage::join_partitions()
{
    /*Identifies the two last partitions of the complete descending path*/
    unsigned int part_one_id = complete_desc_path.at(complete_desc_path.size()-1);
    unsigned int part_two_id = complete_desc_path.at(complete_desc_path.size()-2);
    al_node_t* part_one = partitions.at(part_one_id);
    al_node_t* part_two = partitions.at(part_two_id);

    unsigned int vertex;

    /*Removes the two last partitions of the complete descending path*/
    complete_desc_path.pop_back();
    complete_desc_path.pop_back();

    /*Creates a new average linkage node that combines the partitions of the
    * two partitions being joined*/
    al_node_t* al_node = new al_node_t;
    al_node->partition.reserve(part_one->partition.size()+part_two->partition.size());

    for(unsigned int v = 0; v < part_one->partition.size(); v++)
    {
        vertex = part_one->partition.at(v);
        al_node->partition.push_back(vertex);
    }

    for(unsigned int v = 0; v < part_two->partition.size(); v++)
    {
        vertex = part_two->partition.at(v);
        al_node->partition.push_back(vertex);
    }

    al_node->difference = 0;
    al_node->average = 0;

    /*Sets the average of the new average linkage node*/
    /*
    al_node->average = (double)(part_one->partition.size()*part_one->average +
        part_two->partition.size()*part_two->average) /
        (part_one->partition.size() + part_two->partition.size());
    */
    /*Links the new node to the joined ones*/
    al_node->left = part_one;
    al_node->right = part_two;

    al_node->size = part_one->size + part_two->size;

    /*Don't need to keep the partitions in case they are not leaf ones*/
    if(part_one->partition.size() > 1)
    {
        part_one->partition.clear();
        part_one->partition.reserve(0);
    }

    if(part_two->partition.size() > 1)
    {
        part_two->partition.clear();
        part_two->partition.reserve(0);
    }

    /*The new partition will replace the first joined one, the other will be
     *inactivated*/
    partitions.at(part_one_id) = al_node;
    active_partitions.at(part_two_id) = false;

    /*Updates the distance matrix*/
    double distance;

    for(unsigned int p = 0; p < partitions.size(); p++)
    {
        if(active_partitions.at(p))
        {
            if(p != part_one_id)
            {
                distance = distance_al_nodes_t(partitions.at(p),
                    partitions.at(part_one_id), graph);
            }
            else
            {
                distance = std::numeric_limits<double>::max();
            }

            distance_matrix[p][part_one_id] = distance;
            distance_matrix[part_one_id][p] = distance;
        }
    }
}

/**
 * Prints the average linkage tree recursively
 * @param al_node average linkage node
 * @param depth depth
 * @return
 * @throws
**/
void print_al_node_t(al_node_t* al_node, unsigned int depth)
{
    std::cout << "avg = " << al_node->average << " diff = " << al_node->difference << " size = " << al_node->size << std::endl;

    if(al_node->left == NULL and al_node->right == NULL)
    {
        std::cout << "vertex = " << al_node->partition.at(0) << std::endl;
    }

    if(al_node->left != NULL)
    {
        for(unsigned int d = 0; d <= depth; d++)
        {
            std::cout << " ";
        }

        print_al_node_t(al_node->left, depth+1);
    }

    if(al_node->right != NULL)
    {
        for(unsigned int d = 0; d <= depth; d++)
        {
            std::cout << " ";
        }

        print_al_node_t(al_node->right, depth+1);
    }
}

/**
 * Prints the average linkage tree
 * @param st_node parent node
 * @return
 * @throws
**/
void AverageLinkage::print()
{
    print_al_node_t(tree, 0);
}

/**
 * Frees the memory of the tree recursively
 * @param st_node parent node
 * @return
 * @throws
**/
void free_al_node_t(al_node_t* al_node)
{
    if(al_node->left != NULL)
    {
        free_al_node_t(al_node->left);
    }

    if(al_node->right != NULL)
    {
        free_al_node_t(al_node->right);
    }

    delete al_node;
}

/**
 * Destructor.
 * @param
 * @return
 * @throws
**/
AverageLinkage::~AverageLinkage()
{
    free_al_node_t(tree);
}

/**
 * Writes an average linkage node to a binary file.
 * @param al_node average linkage node
 * @param output_file output file
 * @param node_id variable to define an id for nodes
 * based on the visiting order
 * @return
 * @throws
**/
void write_al_node(al_node_t* al_node, std::ofstream& output_file, unsigned int& node_id)
{
    if(al_node->left != NULL and al_node->right != NULL)
    {
        if(al_node->difference > 0)
        {
            uint32_t id = (uint32_t) node_id;
            output_file.write(reinterpret_cast<char*>(&id), sizeof(uint32_t));
            output_file.write(reinterpret_cast<char*>(&al_node->difference), sizeof(float));
        }

        node_id = node_id + 1;

        write_al_node(al_node->left, output_file, node_id);
        write_al_node(al_node->right, output_file, node_id);
    }
}

/**
 * Writes the wavelet coefficients of the average linkage compression
 * to a file.
 * Format of the binary file:
 * <average dataset><node_id_0><non_zero_difference_0>...
 * <node_id_1><non_zero_difference_1> ...
 * average and differences are (32/64 bits)
 * node_ids are unsigned integers (32 bits) and give the visiting order
 * @param output_file_name ouput file name
 * @return
 * @throws
**/
void AverageLinkage::write(const std::string& output_file_name) const
{
    std::ofstream output_file(output_file_name.c_str(), std::ios::out | std::ios::binary);
    unsigned int node_id = 0;

    /*Writes the overall average to the file*/
    if(tree->average > 0)
    {
        uint32_t id = (uint32_t) node_id;
        output_file.write(reinterpret_cast<char*>(&id), sizeof(uint32_t));
        output_file.write(reinterpret_cast<char*>(&(tree->average)), sizeof(float));
    }

    node_id = node_id + 1;
    write_al_node(tree, output_file, node_id);

    output_file.close();
}

/**
 * Frees the memory of the tree recursively
 * @param wavelets_node_node parent node
 * @return
 * @throws
**/
void free_wavelets_node_t(wavelets_node_t* wavelets_node)
{
    if(wavelets_node->left != NULL)
    {
        free_wavelets_node_t(wavelets_node->left);
    }

    if(wavelets_node->right != NULL)
    {
        free_wavelets_node_t(wavelets_node->right);
    }

    delete wavelets_node;
}

/**
 * Destructor.
 * @param
 * @return
 * @throws
**/
Wavelets::~Wavelets()
{
    free_wavelets_node_t(tree);
}


/**
 * Builds a wavelet tree from the wavelet nodes
 * @param wavelets_nodes wavelets nodes
 * @return
 * @throws
**/
void Wavelets::build_wavelet_tree_recursive(std::vector<wavelets_node_t*>& wavelets_nodes)
{
    wavelets_node_t* new_node;

    if(wavelets_nodes.size() > 1)
    {
        unsigned int node;
        for(node = 0; node < wavelets_nodes.size() - 1; node+=2)
        {
            new_node = new wavelets_node_t;
            new_node->average = 0;
            new_node->difference = 0;
            new_node->size = wavelets_nodes.at(node)->size
                + wavelets_nodes.at(node+1)->size;
            new_node->left = wavelets_nodes.at(node);
            new_node->right = wavelets_nodes.at(node+1);
            new_node->vertex = UINT_MAX;

            wavelets_nodes.at(node/2) = new_node;
        }

        if(node < wavelets_nodes.size())
        {
            wavelets_nodes.at(node/2) = wavelets_nodes.at(node);
            node = node + 2;
        }

        while(wavelets_nodes.size() >= node/2 + 1)
        {
            wavelets_nodes.pop_back();
        }

        build_wavelet_tree_recursive(wavelets_nodes);
    }
}

/**
 * Gets a vector with the intermediate nodes of the wavelet tree
 * recursively
 * @param wavelets_node wavelets node
 * @param wavelets_nodes vector where the nodes will be inserted
 * @return
 * @throws
**/
void get_wavelets_intermediate_nodes(wavelets_node_t* wavelets_node,
    std::vector<wavelets_node_t*>& wavelets_nodes)
{
    if(wavelets_node->left != NULL && wavelets_node->right != NULL)
    {
        wavelets_nodes.push_back(wavelets_node);
        wavelets_node->left->average = 0;
        wavelets_node->right->average = 0;
        get_wavelets_intermediate_nodes(wavelets_node->left, wavelets_nodes);
        get_wavelets_intermediate_nodes(wavelets_node->right, wavelets_nodes);
    }
}

/**
 * Wavelet coefficient pruning
 * @param
 * @return
 * @throws
**/
void Wavelets::keep_top_coefficients()
{
    std::vector<wavelets_node_t*> wavelets_nodes;
    wavelets_nodes.reserve(graph->size());
    unsigned int num_remaining_coefficients = num_coefficients;

    /*Gets the vector of wavelets nodes*/
    get_wavelets_intermediate_nodes(tree, wavelets_nodes);

    /*Sorts the wavelets nodes*/
    std::sort(wavelets_nodes.begin(), wavelets_nodes.end(), CompareWaveletsNodes());

    /*Decides whether the global average, which is a kind of special
     * coefficient will be kept, the relevance of the average coefficient*/
    if(num_coefficients == 0 ||
        wavelets_nodes.at(num_coefficients-1)->difference*2
        > tree->average * tree->size)
    {
        tree->average = 0;
        num_remaining_coefficients = num_coefficients;
    }
    else
    {
        num_remaining_coefficients = num_remaining_coefficients - 1;
    }

    /*Sets the coefficients not kept to 0*/
    for(unsigned int c = num_remaining_coefficients;
        c < wavelets_nodes.size(); c++)
    {
        wavelets_nodes.at(c)->difference = 0;
    }
}

Wavelets::Wavelets(Graph& _graph):
    GraphCompressionAlgorithm(_graph)
{;}

void Wavelets::build_wavelet_tree()
{
    std::vector<wavelets_node_t*> wavelets_nodes;
    wavelets_nodes.reserve(graph->size());

    wavelets_node_t* wavelets_node;

    for(unsigned int v = 0; v < graph->size(); v++)
    {
        wavelets_node = new wavelets_node_t;
        wavelets_node->average = 0;
        wavelets_node->difference = 0;
        wavelets_node->size = 1;
        wavelets_node->vertex = graph->at(v);
        wavelets_node->left = NULL;
        wavelets_node->right = NULL;

        wavelets_nodes.push_back(wavelets_node);
    }

    build_wavelet_tree_recursive(wavelets_nodes);

    tree = wavelets_nodes.at(0);
}

void compute_average_coefficients_wavelets_node_t(wavelets_node_t* wavelets_node, Graph* graph)
{
    if(wavelets_node->left != NULL and wavelets_node->right != NULL)
    {
        compute_average_coefficients_wavelets_node_t(wavelets_node->left, graph);
        compute_average_coefficients_wavelets_node_t(wavelets_node->right, graph);

        wavelets_node->average =
            (double)(wavelets_node->left->size
            * wavelets_node->left->average
            + wavelets_node->right->average
            * wavelets_node->right->size)
            / (wavelets_node->left->size
            + wavelets_node->right->size);
    }
    else
    {
        wavelets_node->average = graph->value(wavelets_node->vertex);
    }
}

void compute_difference_coefficients_wavelets_node_t(wavelets_node_t* wavelets_node, Graph* graph)
{
    if(wavelets_node->left != NULL and wavelets_node->right != NULL)
    {
        compute_difference_coefficients_wavelets_node_t(wavelets_node->left, graph);
        compute_difference_coefficients_wavelets_node_t(wavelets_node->right, graph);

        wavelets_node->difference =
            (double) (wavelets_node->average
            - wavelets_node->left->average)
            * wavelets_node->left->size;
    }
}

/**
 * Computes the average coefficients for a wavelet
 * tree
 * @param
 * @return
 * @throws
**/
void Wavelets::compute_average_coefficients()
{
    compute_average_coefficients_wavelets_node_t(tree, graph);
}

/**
 * Computes the difference coefficients for a wavelet
 * tree
 * @param
 * @return
 * @throws
**/
void Wavelets::compute_difference_coefficients()
{
    compute_difference_coefficients_wavelets_node_t(tree, graph);
}

/**
 * Runs the slice tree compression.
 * @param
 * @return
 * @throws
**/
void Wavelets::compress(const unsigned int budget)
{
    budget_compression = budget;
    num_coefficients = (unsigned int) floor(budget_compression
        / SIZE_FLOAT_INT);

    build_wavelet_tree();
    compute_average_coefficients();
    compute_difference_coefficients();
    keep_top_coefficients();
}


const double Wavelets::compute_sse()
{
    double sse_value = 0;

    for(unsigned int v = 0; v < graph->size(); v++)
    {
        sse_value = sse_value + pow(graph->orig_value(v) - value(v), 2     );
    }

    return sse_value;
}

/**
 * Prints a wavelet tree recursively.
 * @param wavelets_node wavelets tree node
 * @param depth depth
 * @return
 * @throws
**/
void print_wavelets_node_t(wavelets_node_t* wavelets_node, unsigned int depth)
{
    std::cout << "avg = " << wavelets_node->average << " diff = " << wavelets_node->difference << " size = " << wavelets_node->size << std::endl;

    if(wavelets_node->left != NULL)
    {
        for(unsigned int d = 0; d <= depth; d++)
        {
            std::cout << " ";
        }

        print_wavelets_node_t(wavelets_node->left, depth+1);
    }

    if(wavelets_node->right != NULL)
    {
        for(unsigned int d = 0; d <= depth; d++)
        {
            std::cout << " ";
        }

        print_wavelets_node_t(wavelets_node->right, depth+1);
    }
}

/**
 * Prints a wavelet tree.
 * @param
 * @return
 * @throws
**/
void Wavelets::print()
{
    print_wavelets_node_t(tree, 0);
}

/**
 * Writes a wavelets node to a file in preorder
 * @param wavelets_node wavelets tree node
 * @param output_file output file
 * @param node_id id of the coefficient
 * @return
 * @throws
**/
void write_wavelets_node(wavelets_node_t* wavelets_node, std::ofstream& output_file,
    unsigned int& node_id)
{
    if(wavelets_node->left != NULL && wavelets_node->right != NULL)
    {
        if(wavelets_node->difference > 0)
        {
            uint32_t id = (uint32_t) node_id;
            output_file.write(reinterpret_cast<char*>(&id), sizeof(uint32_t));
            output_file.write(reinterpret_cast<char*>(&wavelets_node->difference),
                sizeof(float));
        }

        node_id = node_id + 1;

        write_wavelets_node(wavelets_node->left, output_file, node_id);
        write_wavelets_node(wavelets_node->right, output_file, node_id);
    }
}

/**
 * Writes the wavelet coefficients to a binary file.
 * Format of the binary file:
 * <average dataset><non_zero_coefficient_id_1>
 * <non_zero_coefficient_value_1> ...
 * difference coefficients are floats (32/64 bits)
 * coefficient ids are unsigned ints (32/64 bits)
 * @param output_file_name ouput file name
 * @return
 * @throws
**/
void Wavelets::write(const std::string& output_file_name) const
{
    std::ofstream output_file(output_file_name.c_str(), std::ios::out | std::ios::binary);
    unsigned int node_id = 0;

    /*Writes the overall average to the file*/
    if(tree->average > 0)
    {
        uint32_t id = (uint32_t) node_id;
        output_file.write(reinterpret_cast<char*>(&id), sizeof(uint32_t));
        output_file.write(reinterpret_cast<char*>(&(tree->average)), sizeof(float));
    }

    node_id = node_id + 1;
    write_wavelets_node(tree, output_file, node_id);

    output_file.close();
}

/**
 * Sets the values (recursively) as they would be recovered after compression
 * using wavelets.
 * @param wavelets_node wavelets tree node
 * @param values vector in which the values will be set
 * @return
 * @throws
**/
void set_compressed_values_wavelets_node_t(wavelets_node_t* wavelets_node,
    std::vector<double>& values)
{
    if(wavelets_node->left != NULL and wavelets_node->right != NULL)
    {
        wavelets_node->left->average = wavelets_node->average
            - (double) wavelets_node->difference / wavelets_node->left->size;

        wavelets_node->right->average = wavelets_node->average
            + (double) wavelets_node->difference / wavelets_node->right->size;

        set_compressed_values_wavelets_node_t(wavelets_node->left, values);
        set_compressed_values_wavelets_node_t(wavelets_node->right, values);
    }
    else
    {
        values.at(wavelets_node->vertex) = wavelets_node->average;
    }
}

/**
 * Sets the values as they would be recovered after compression
 * using wavelets.
 * @param
 * @return
 * @throws
**/
void Wavelets::set_compressed_values()
{
    values.clear();

    for(unsigned int v = 0; v < graph->size(); v++)
    {
        values.push_back(0);
    }

    set_compressed_values_wavelets_node_t(tree, values);
}

