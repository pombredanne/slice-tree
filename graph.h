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
 * Definitions of a class for graph manipulation
**/

#ifndef GRAPH_H
#define GRAPH_H

#define MAX_GRAPH_DIAMETER UCHAR_MAX
#define MAX_SIZE_DISTANCE_STR 2000000000

/*std includes*/
#include <string>
#include <exception>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <map>
#include <algorithm> 
#include <sstream>
#include <cmath>
#include <pthread.h>

typedef struct Edge
{
	unsigned int v_one;
	unsigned int v_two;
	double difference;
}edge_t;

class ComparePairs
{
	public:
		bool operator()(const std::pair<unsigned int, double>* p_one, 
			const std::pair<unsigned int, double>* p_two) const
		{
			return p_one->second < p_two->second;
		}
};

typedef struct _PthreadParameters
{
	std::vector< std::list<unsigned int>* >* adjacency_list;
	std::list<unsigned int>* pool;
	pthread_mutex_t* mutex_pool;
	std::vector<std::vector<unsigned int>*>* partition_sizes;
	unsigned int max_radius;
}PthreadParameters;

/**
 * Class for graph manipulation
**/
class Graph
{
	public:
		/**
		 * Graph constructor
		 * @param  graph_file_name input graph file name
		 * @param  value_file_name input values file name
		 * @return 
		 * @throws 
		 **/

		Graph(const std::string& graph_file_name, const std::string& values_file_name, const bool directed) throw (std::ios_base::failure);
		
		/**
		 * Graph constructor
		 * @param  graph_file_name input graph file name
		 * @return 
		 * @throws 
		 **/
		Graph(const std::string& graph_file_name) throw (std::ios_base::failure);
		
		/**
		 * Graph destructor
		 * @param 
		 * @return 
		 * @throws 
		**/
		virtual ~Graph();
		
		/**
		 * Builds a distance structure for slice tree, this structure 
		 * efficiently returns the list of vertices at a given distance
		 * from a certain vertex
		 * @param 
		 * @return 
		 * @throws 
		**/
		void build_distance_str_slice_tree();
		
		/**
		 * Builds a distance structure for slice tree using a set of sample 
		 * vertices, this strcutre efficiently returns the list of vertices 
		 * at a given distance from a certain vertex
		 * @param max_radius maximum radius for slice tree
		 * @return 
		 * @throws 
		**/
		void build_distance_str_slice_tree_sample
			(const unsigned int max_radius);
		
		void build_distance_str_slice_tree_sample
			(const unsigned int max_radius, 
			std::vector<unsigned int>& partition);
		
		void start_distance_str_slice_tree_sample();

		void clear_distance_str_sample(const std::vector<unsigned int>& partition,
			const std::vector<bool>& bitmap);

		void clear_distance_str_biased_sample(const unsigned int num_samples,
			const std::vector<unsigned int>& partition);
		
		/**
		 * Prints the slice tree distance structure
		 * @param 
		 * @return 
		 * @throws 
		**/
		void print_distance_str_slice_tree();

		/**
		 * Builds a distance matrix for the graph
		 * @param 
		 * @return 
		 * @throws 
		**/
		void build_distance_matrix();

		/**
		 * Builds a BFS-based sorted vector for the graph
		 * @param 
		 * @return 
		 * @throws 
		**/
		void build_bfs_vector();

		/**
		 * Builds a priority first-based sorted vector for the graph
		 * @param 
		 * @return 
		 * @throws 
		**/
		void build_priority_first_vector(const unsigned int start);
		
		/**
		 * Performs a BFS search from the center and bounded by radius distance
		 * including in visited only the vertices in the bitmap
		 * @param visited vector updated with the vertices visited in the search
		 * @param center center
		 * @param radius radius
		 * @param bitmap bitmap for the partition
		 * @return 
		 * @throws 
		**/
		void bounded_bfs(std::vector<unsigned int>& visited, const unsigned int center, const unsigned int radius, const std::vector<bool>& bitmap) const;

		/**
		 * Updates some data structures used to compute upper bounds on 
		 * the sizes of partitions generated by a slice.
		 * @param center center
		 * @param radius radius
		 * @param dist_near_center distance to the center of the nearest
		 * partition for each each vertex.
		 * @param dist_center_part distance to the center of the partition
		 * where the slice lies in for each vertex.
		 * @param radius_near_center radius of the nearest partition
		 * for each vertex
		 * @param radius_part radius of the partition the slice lies in.
		 * @throws
		 * @return
		**/
		void update_partition_size_structs(
			const unsigned int center,
			const unsigned int radius,
			std::vector<unsigned int>& dist_near_center,
			std::vector<unsigned int>& dist_center_part,
			std::vector<unsigned int>& radius_near_center,
			std::vector<unsigned int>& radius_part,
			const std::vector<bool>& bitmap) const;

		/**
		 * Selects a set of vertices as random samples 
		 * (with replacement) from the graph
		 * @param num_samples number of samples
		 * @param partition partition to be sampled from 
		 * @return
		 * @throws 
		**/
		void set_uniform_sample(const unsigned int num_samples,
			const std::vector<unsigned int>& partition);

		/**
		 * Selects a set of vertices as random samples 
		 * (with replacement) from the graph
		 * @param num_samples number of samples
		 * @param partition partition to be sampled from 
		 * @return
		 * @throws 
		**/
		void resample_uniform_sample(const unsigned int num_samples,
			const std::vector<unsigned int>& partition);

		/**
		 * Selects (with replacement) a set of vertices from the graph 
		 * in a biased way  where the bias is proportional to 
		 * |value(vertex)-mean(partition)|.
		 * @param num_samples number of samples
		 * @param partition partition to be sampled from 
		 * @return
		 * @throws
		**/
		void set_biased_sample(const unsigned int num_samples, 
			const std::vector<unsigned int>& partition);

		/**
		 * Selects (with replacement) a set of vertices from the graph 
		 * in a biased way  where the bias is proportional to 
		 * |value(vertex)-mean(partition)|.
		 * @param num_samples number of samples
		 * @param partition partition to be sampled from 
		 * @return
		 * @throws
		**/
		void resample_biased_sample(const unsigned int _num_samples,
			const std::vector<unsigned int>& partition);
		
		/**
		 * Selects a sample from the graph (biased or unbiased).
		 * @param num_samples number of samples
		 * @param partition partition to be sampled from 
		 * @return
		 * @throws
		**/
		void set_sample(const unsigned int num_samples,
			const std::vector<unsigned int>& partition);
		
		/**
		 * Selects a sample from the graph (biased or unbiased).
		 * @param num_samples number of samples
		 * @param partition partition to be sampled from 
		 * @return
		 * @throws
		**/
		void resample(const unsigned int num_samples,
			const std::vector<unsigned int>& partition);
		/**
		 * Prints the graph (for debugging purposes)
		 * @param
		 * @return
		 * @throws 
		**/
		void print();

		/**
		 * Builds the distance structure for a single center
		 * @param center center
		 * @param vertices_at_distance structure
		 * @param radius max radius to be searched
		 * @return
		 * @throws
		**/
		void build_distance_str_slice_tree_vertex(unsigned int center, 
			std::vector<std::list<unsigned int>*>& 
			vertices_at_distance, 
			const unsigned int max_radius) const;
		
		/**
		 * Reads the pre-computed partition sizes from a file.
		 * @param input_file_name input file
		 * @param max_radius maximum radius
		 * @return
		 * @throws
		**/
		void read_partition_sizes(const std::string& input_file_name,
			const unsigned int max_radius);
		
		/**
		 * Pre-computes the partitions sizes for all centers and radius
		 * and saves them in a file. Uses multiple threads to speed up
		 * the computations.
		 * @param num_threads number of threads available
		 * @param output_file_name output file
		 * @param max_radius maximum radius
		 * @return
		 * @throws
		**/
		void pre_compute_partition_sizes(const unsigned int num_threads, 
			const std::string& output_file_name,
			const unsigned int max_radius);
		
		/*Inline methods:*/
		/**
		 * Returns the size of the graph
		 * @param 
		 * @return size
		 * @throws 
		**/
		const inline unsigned int size() const
		{
			return adjacency_list.size();
		}

		/**
		 * Returns the diameter of the graph
		 * @param 
		 * @return diameter
		 * @throws 
		**/
		const inline unsigned int diameter() const
		{
			return graph_diameter;
		}

		/**
		 * Returns the number of vertices at a given distance 
		 * from a vertex given as parameter
		 * @param vertex vertex
		 * @param distance distance
		 * @return number of vertices
		 * @throws 
		**/
		inline std::list<unsigned int>* vertices_at_distance(const unsigned int vertex, 
			const unsigned int distance) const
		{
//			printf("distance = %d, size = %d, vertex = %d, size = %d\n", 
//				distance, distance_str[vertex]->size(), vertex, distance_str.size());
			return distance_str[vertex]->at(distance);
		}

		/**
		 * Returns the distance between two vertices
		 * @param v_one vertex
		 * @param v_two vertex
		 * @return distance
		 * @throws 
		**/
		const inline unsigned int distance(const unsigned int v_one, 
			const unsigned int v_two) const
		{
			return distance_matrix[v_one][v_two];
		}

		/**
		 * Returns the maximum distance from a given vertex 
		 * to any other vertex
		 * @param vertex
		 * @return maximum distance
		 * @throws 
		**/
		const inline unsigned int max_distance(const unsigned int vertex) const 
		{
			return distance_str[vertex]->size();
		}
		
		/**
		 * Returns the value of a vertex
		 * @param 
		 * @return value
		 * @throws 
		**/
		const inline double value(const unsigned int v) const
		{
			if(biased_sampling)
			{
				return (double) (vertex_values.at(v) * lambda) / fabs(vertex_values.at(v) - mu);
			}
			else
			{
				return vertex_values.at(v);
			}
		}

		/**
		 * Returns the value of lambda
		 * @param 
		 * @return lambda
		 * @throws 
		**/
		const inline double get_lambda()
		{
			return lambda;
		}
		
		/**
		 * Returns the weight of a vertex
		 * @param v vertex
		 * @return weight
		 * @throws 
		**/
		const inline double weight(const unsigned int v) const
		{
			if(biased_sampling)
			{
				return (double) lambda / fabs(vertex_values.at(v) - mu);
			}
			else
			{
				return 1;
			}
		}
		
		/**
		 * Returns the actual value of a vertex
		 * @param v vertex
		 * @return weight
		 * @throws 
		**/
		const inline double orig_value(const unsigned int v) const
		{
			return vertex_values.at(v);
		}

		/**
		 * Returns the name of a vertex
		 * @param 
		 * @return diameter
		 * @throws 
		**/
		const inline std::string name(const unsigned int v) const
		{
			return vertex_names.at(v);
		}

		/**
		 * Sets a value to a vertex in the graph
		 * @param v vertex
		 * @param value value
		 * @return 
		 * @throws 
		**/
		inline void set_value(const unsigned int v, const double value)
		{
			vertex_values[v] = value;
		}
		
		/**
		 * Returns the id of the vertex at the pos position 
		 * of the sorted vector
		 * @param pos position
		 * @return vertex id
		 * @throws 
		**/
		const inline unsigned int at(const unsigned int pos) const
		{
			return sorted_vector.at(pos);
		}

		/**
		 * Returns the number of occurrences of a given item in the sample 
		 * @param vertex id
		 * @return count
		 * @throws 
		**/
		const inline unsigned int count(const unsigned int v) const
		{
			if(biased_sampling or uniform_sampling) return count_sample[v]; else return 1;
		}

		/**
		 * Sets the sampling to biased, but the samples are not
		 * produced again.
		 * @param 
		 * @return
		 * @throws 
		**/
		void inline set_biased_sampling()
		{
			biased_sampling = true;
			uniform_sampling = false;
		}
		
		/**
		 * Sets the sampling to uniform, but the samples are not
		 * produced again.
		 * @param 
		 * @return
		 * @throws 
		**/
		void inline set_uniform_sampling()
		{
			uniform_sampling = true;
			biased_sampling = false;
		}

		/**
		 * Returns the partition size for a given center and radius
		 * @param center center 
		 * @param radius radius
		 * @return partition size
		 * @throws 
		**/
		const inline unsigned int get_partition_size(const unsigned int center, 
			const unsigned int radius) const
		{
			if(radius < partition_sizes.at(center)->size())
			{
				return partition_sizes.at(center)->at(radius);
			}
			else
			{
				return partition_sizes.at(center)->back();
			}
		}

		/**
		 * Returns the sum of values
		 * @param
		 * @return sum of values
		 * @throws
		**/
		const inline double get_sum_values()
		{
			return sum_values;
		}

		/**
		 * Returns the sum of weights
		 * @param
		 * @return sum of weights
		 * @throws
		**/
		const inline double get_sum_weights()
		{
			return sum_weights;
		}

		const inline double get_sum_weighted_values()
		{
			return sum_weighted_values;
		}

		const inline double get_theta()
		{
			return theta;
		}
	private:
		std::vector< std::list<unsigned int>* > adjacency_list;
		std::vector< std::list<unsigned int>* > back_adjacency_list;
		unsigned short int** distance_matrix;
		std::vector<double> vertex_values;
		std::map<unsigned int, std::string> vertex_names;
		std::map<std::string,unsigned int> vertex_ids;
		std::vector< std::vector< std::list<unsigned int >* >* > distance_str;
		std::vector<unsigned int> sorted_vector; 
		unsigned int graph_diameter;
		std::vector<unsigned int> count_sample;
		std::list<unsigned int> samples;
		double lambda;
		double mu;
		bool biased_sampling;
		bool uniform_sampling;
		double sum_values;
		double sum_weights;
		double sum_weighted_values;
		bool directed;
		std::vector<std::vector<unsigned int>*> partition_sizes;
		std::vector<bool> is_sampled;
		double theta;
		std::vector<float> selection_prob;
		std::vector<unsigned int> distances;
		unsigned long size_distance_str;
		
		/**
		 * Performs a bfs search over the graph, returning the size of the set of vertices 
		 * visited and updating the visited set.
		 * @param root root
		 * @param visited visited set
		 * @return 
		 * @throws 
		**/
		const unsigned int bfs(const unsigned root, std::vector<bool>& visited);
		
		/**
		 * Reads the graph data
		 * @param graph_file_name input file with edges
		 * @param values_file_name input file with vertex values
		 * @return
		 * @throws ios_base::failure in case can't read the input files
		 **/
		void read_graph(const std::string& graph_file_name, const std::string& values_file_name) throw (std::ios_base::failure);

		/**
		 * Reads the graph data
		 * @param graph_file_name input file with edges
		 * @return
		 * @throws ios_base::failure in case can't read the input files
		 */
		void read_graph(const std::string& graph_file_name) throw (std::ios_base::failure);

		/**
		 * Counts the number of vertices in the input graph
		 * @param values_file_name file with the vertex values
		 * @return number of vertices
		 * @throws
		**/
		unsigned int count_vertices(const std::string& values_file_name) const;

		void free_distance_str();
};

#endif
