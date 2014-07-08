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
 * Definitions of a class for graph compression
**/

#ifndef GRAPHCOMPRESSION_H
#define GRAPHCOMPRESSION_H

/*std includes*/
#include <string>
#include <exception>
#include <list>
#include <vector>
#include <cmath>
#include <algorithm>
#include <climits>
#include <limits>
#include <sstream>

/*my includes*/
#include "graph.h"
#include "perf.h"

#define SIZE_FLOAT_INT 8

/**
 * Slice tree node
**/ 
typedef struct STNode
{
	/*
	 * Cut performed in case the partition was split
	 * In case the partition was not split, stores
	 * the best candidate split information
	**/
	unsigned int center;
	unsigned int radius;
	
	unsigned int diameter;
	unsigned int size;
	
	float average;
	float difference;

	/*children*/
	struct STNode* left;
	struct STNode* right;
	
	/*Error of the partition at this level of the tree*/
	double error_partition;
	
	/*Error of the best cut i.e. for the center/radius defined*/
	double error_best_cut;
	
	/*Vertices in the partition*/
	std::vector<unsigned int> partition;

	/** Bitmap for efficiently checking whether a vertex is part
	 * of the partition.
	**/
	std::vector<bool> in_partition;
}st_node_t;

/**
 * Probabilistic upper-bound
**/
typedef struct UpperBoundType
{
	unsigned int center;
	unsigned int radius;
	double bound;
	double estimate;
	std::list< std::pair<unsigned int, double>* > bounds;
}up_bound_t;

/**
 * Prints a slice tree node
 * @param st_node slice tree node
 * @param depth depth of the node
 * @return 
 * @throws
**/
void print_st_node(st_node_t* st_node, unsigned int depth, std::string pid, std::string type, Graph* graph);

/**
 * Compares two slice tree nodes
 * By using this function you get an increasing order.
**/
class CompareCuts
{
	public:
		bool operator()(const st_node_t* n_one, const st_node_t* n_two) const
		{
			return n_one->error_best_cut > n_two->error_best_cut;
		}
};

/**
 * Generic class that implements a compression algorithm
**/
class GraphCompressionAlgorithm
{
	public:
		/**
		 * Constructor.
		 * @param graph_to_compress input graph
		 * @param budget available budget
		 * @return 
		 * @throws
		**/
		GraphCompressionAlgorithm(Graph& graph_to_compress);
		
		/**
		 * Constructor.
		 * @param input_file_name input file with a serialized compressed graph
		 * @param graph graph for which the values will be decompressed
		 * @return 
		 * @throws
		**/
		GraphCompressionAlgorithm(const std::string& input_file_name, Graph& graph);
		
		/**
		 * Destructor. Does nothing. Should be defined by the child class.
		 * @param
		 * @return 
		 * @throws
		**/
		virtual ~GraphCompressionAlgorithm(){};
		
		/**
		 * Does nothing. Should be defined by the child class.
		 * @param
		 * @return 
		 * @throws
		**/
		virtual void compress(const unsigned int budget){};

		/**
		 * Does nothing. Should be defined by the child class.
		 * @param
		 * @return 
		 * @throws
		**/
		virtual void decompress(){};

		/**
		 * Does nothing. Should be defined by the child class.
		 * @param
		 * @return 
		 * @throws
		**/
		virtual void write(const std::string& output_file_name) const{};

		/**
		 * Does nothing. Should be defined by the child class.
		 * @param 
		 * @return
		 * @throws
		**/
		virtual void set_compressed_values(){};
		
		/**
		 * Returns the value of a vertex as it would be compressed.
		 * @param
		 * @return 
		 * @throws
		**/
		inline const double value(unsigned int vertex)
		{
			return values[vertex];
		}

		/**
		 * Sets values to the graph based on recovered slice tree
		 * @param 
		 * @return 
		 * @throws
		**/
	protected:
		Graph* graph;
		unsigned int budget_compression;
		std::string compressed_file_name;
		std::vector<double> values;
};

/**
 * Class that implements the slice tree compression
**/
class SliceTree: public GraphCompressionAlgorithm
{
	public:
		/**
		 * Constructor. 
		 * @param graph graph 
		 * @param max_radius maximum radius for slice tree
		 * @param exhaustive_split consider all splits in all partitions if set
		 * @return 
		 * @throws 
		**/
		SliceTree(Graph& graph, const unsigned int _max_radius, const bool _exhaustive_split):
			GraphCompressionAlgorithm(graph)
		{
			max_radius = _max_radius;
			exhaustive_split = _exhaustive_split;
		}

		/**
		 * Constructor. 
		 * @param input_file_name serialized graph with compressed data 
		 * @param graph graph
		 * @return 
		 * @throws 
		**/
		SliceTree(const std::string& input_file_name, Graph& graph):
			GraphCompressionAlgorithm(input_file_name, graph){;}

		/**
		 * Builds a slice tree from the serialized content of a file
		 * @param input_file_name input file
		 * @return
		 * @throws
		**/
		void decompress();

		/**
		 * Wrapper for compress heuristics
		 * @param
		 * @return 
		 * @throws
		**/
		void compress(const unsigned int budget);
		
		/**
		 * Runs the slice tree compression. All partitions probed
		 * @param
		 * @return 
		 * @throws
		**/
		void compressExhaustive(const unsigned int budget); 

		/**
		 * Runs the slice tree compression. Only the highest error slice is
 		 * probed for all possible cuts 
		 * @param
		 * @return 
		 * @throws
		**/
		void compressGreedy(const unsigned int budget);

		/**
		 * Destructor
		 * @param
		 * @return
		 * @throws
		**/
		virtual ~SliceTree();

		/**
		 * Writes the slice tree to a file
		 * @param output_file_name ouput file name
		 * @return
		 * @throws
		**/
		void write(const std::string& output_file_name) const;

		/**
		 * Prints the slice tree on the terminal
		 * @param 
		 * @return
		 * @throws
		**/
		void print() const;
		
		/**
		 * Returns the budget of the compression
		 * @param num_partitions number of partitions
		 * @param num_vertices number of vertices in the graph
		 * @param diameter diameter of the graph
		 * @return number of partitions
		 * @throws
		 **/
		const static unsigned int budget(const unsigned int num_partitions, Graph& graph)
		{
			return SIZE_FLOAT_INT + (num_partitions - 1) * 
				size_node(graph.size(), graph.diameter());
		}

		/**
		 * Sets the values as they would be recovered after compression
		 * using slice tree.
		 * @param 
		 * @return
		 * @throws
		**/
		void set_compressed_values();

	protected:
		st_node_t* tree;
		unsigned int n_partitions;
		double global_error; //Keeps the final sse
		unsigned int max_radius;
		bool exhaustive_split;

		/**
		 * Extends the slice tree recovered from a serialized file
		 * @param 
		 * @return
		 * @throws
		**/
		void extend_tree();
		void extend_st_node(st_node_t* st_node, Graph* graph);

		/**
		 * Computes the sse of a partition
		 * @param partition partition
		 * @return sse
		 * @throws
		**/
		const double sse_partition(const std::vector<unsigned int>& partition) const;
		
		/**
		 * Find partition of max SSE
		 * @param
		 * @return 
		 * @throws
		**/
		st_node_t* getMaxSSEPartiton(st_node_t* root);
		
		/**
		 * Identifies the optimal cut (center/radius) for the
		 * partition represented as a slice tree node
		 * @param st_node slice tree node
		 * @return 
		 * @throws
		**/
		virtual void optimal_cut(st_node_t* st_node);
		
		/**
		 * Identifies the optimal radius for a given center and partition
		 * @param center center to be considered
		 * @partition partition to be split
		 * @diameter diameter of the partition to be split
		 * @parameter in_partition bitmap of the partition to be split
		 * @parameter average average value of the partition
		 * @return pair <error, radius>
		 * @throws
		**/
		const std::pair<double, unsigned int> min_error_radius(const unsigned int center, const std::vector<unsigned int>& partition, unsigned int diameter, const std::vector<bool>& in_partition, const double average) const;
		
		/**
		 * Computes the average value of a partition
		 * @param partition partition
		 * @return average
		 * @throws
		**/
		const double average_partition(const std::vector<unsigned int>& partition) const;
		
		/**
		 * Splits a partition given the center and radius defined
		 * by the slice tree node
		 * @param st_node slice tree node
		 * @return true in case the split was performed
		 *	false, otherwise
		 * @throws
		**/
		virtual const bool split_partition(st_node_t* st_node);

		/**
		  * Computes the difference coefficients for a wavelet-like 
		  * decomposition of the slice tree
		  * @param
		  * @return
		  * @throws
		 **/
		void compute_difference_coefficients();
		
		/**
		 * Clears all the partition in the slice tree 
		 * recursivelly.
		 * @param 
		 * @return
		 * @throws
		**/
		void clear_partitions();
	
		/**
		 * Computes the size of a slice tree node in bytes
		 * It is important to notice that his size is theoretical
		 * in the sense that I'm assuming there is a minimal representation
		 * of a node that is not implemented here.
		 * @param num_vertices number of vertices in the graph
		 * @param diameter diameter of the graph
		 * @return size of the slice tree node
		 * @throws
		 **/
		const static inline unsigned int size_node(const unsigned int num_vertices, 
			const unsigned int diameter)
		{
			return (int) ceil((float)(ceil(log2(num_vertices)) + ceil(log2(diameter+1)) + 8*SIZE_FLOAT_INT + 2) / 8);
		}
		
		/**
		 * Computes the number of partitions of the compression
		 * @param budget budget
		 * @param num_vertices number of vertices in the graph
		 * @param diameter diameter of the graph
		 * @return number of partitions
		 * @throws
		 **/
		const static inline unsigned int num_partitions(const unsigned int budget, 
			const unsigned int num_vertices, const unsigned int diameter)
		{
			return 2 + (int)floor((float) (budget - size_node(num_vertices, diameter) - SIZE_FLOAT_INT) / size_node(num_vertices, diameter));
		}
};

/**
 * Class that implements the slice tree compression using sampling
 * to identify probabilistic good slices.
**/
class SliceTreeSamp: public SliceTree
{
	public:	
		/**
		 * Constructor. 
		 * @param graph graph 
		 * @param max_radius maximum radius for slice tree
		 * @param exhaustive_split consider all splits in all partitions if set		
		 * @param delta probability for bounds in error reduction estimates
		 * for slices
		 * @param sampling_rate sampling rate
		 * @return 
		 * @throws 
		**/
		SliceTreeSamp(Graph& _graph, const unsigned int max_radius, 
			const bool _exhaustive_split, 
			const double _delta, 
			const double _sampling_rate,
			const double _rho):
			SliceTree(_graph, max_radius, _exhaustive_split)
		{
			delta = _delta;
			theta = compute_theta();
			sampling_rate = _sampling_rate;
			rho = _rho;
			_graph.start_distance_str_slice_tree_sample();
			
			/*Initializing data structures that are used to compute
			* upper and lower bounds on the sizes of partitions*/
			dist_near_center.reserve(graph->size());
			radius_near_center.reserve(graph->size());
			dist_center_part.reserve(graph->size());
			radius_part.reserve(graph->size());

			for(unsigned int v = 0; v < graph->size(); v++)
			{
				dist_near_center.push_back(UINT_MAX);
				radius_near_center.push_back(UINT_MAX);
				dist_center_part.push_back(UINT_MAX);
				radius_part.push_back(UINT_MAX);
			}
		}
		
		inline static unsigned int count_bound_one()
		{
			return num_pruned_bound_1;
		}
		
		inline static unsigned int count_bound_two()
		{
			return num_pruned_bound_2;
		}

		inline static unsigned int count_bound_three()
		{
			return num_pruned_bound_3;
		}
		
		inline static double pruned()
		{
			if(total_slices > 0)
			{
				return (double) num_pruned / total_slices;
			}
			else
			{
				return 0;
			}
		}

		/**
		 * Destructor
		 * @param
		 * @return
		 * @throws
		**/
		virtual ~SliceTreeSamp(){;}
	protected:
		double delta;
		double theta;
		double sampling_rate;
		double rho;
		static unsigned int num_pruned_bound_1;
		static unsigned int num_pruned_bound_2;
		static unsigned int num_pruned_bound_3;
		static unsigned int total_slices;
		static unsigned int num_pruned;

		/*Data structures for computing upper and lower
		 * bounds on sizes of partitions without actually
		 * going through the complete list of vertices in 
		 * the slice*/
		std::vector<unsigned int> dist_near_center;
		std::vector<unsigned int> dist_center_part;
		std::vector<unsigned int> radius_near_center;
		std::vector<unsigned int> radius_part;
		
		/**
		 * Identifies a probabilistic optimal cut 
		 * (center/radius) for the partition represented 
		 * as a slice tree node using sampling. This overwrites
		 * the standard function, which makes exact computations.
		 * @param st_node slice tree node
		 * @return
		 * @throws
		**/
		void optimal_cut(st_node_t* st_node);
		void optimal_cut_exact(st_node_t* st_node) const;
		
		/**
		 * Computes theta, which is the range in which all values are.
		 * @param 
		 * @return
		 * @throws
		**/
		double compute_theta();
		
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
		void upper_bound_error_reduction(
			up_bound_t* up_bound,
			const unsigned int center,
		  	const std::vector<unsigned int>& partition,
		        const unsigned int diameter, 
		        const std::vector<bool>& in_partition, 
		        const double average,
			const double sum_weighted_values,
			const double sum_weights,
			const unsigned int total_samples,
			const double sse_partition) const;
		
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
		const unsigned lower_bound_size_partition(
			const unsigned int center,
		        const unsigned int radius,
		        const std::vector<unsigned int>& partition) const;
		
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
		const unsigned int upper_bound_size_partition(const unsigned int center,
		         const unsigned int radius, 
			 const std::vector<unsigned int>& partition) const;
		
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
		const unsigned int lower_bound_size_comp_partition(
			const unsigned int center,
			const unsigned int radius, 
			const std::vector<unsigned int>& partition) const;

		const unsigned int upper_bound_size_comp_partition(
			const unsigned int center,
			const unsigned int radius,
			const std::vector<unsigned int>& partition) const;

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
		std::pair<double, double>
		upper_bound_error_reduction_mean_estimate_out(
			const unsigned int center, const unsigned int radius,
			const std::vector<unsigned int>& partition,
			const double average, const double weighted_mean,
			const unsigned int num_samples_part) const;
		
		/**
		 * Splits a partition given the center and radius defined
		 * by the slice tree node when sampling is applied. The 
		 * difference from the standard version is that here we
		 * update some data structures for computing upper and lower
		 * bounds for sizes of partitions.
		 * @param st_node slice tree node
		 * @return true in case the split was performed
		 *      false, otherwise
		 * @throws
		**/
		virtual const bool split_partition(st_node_t* st_node);
		
		 /**
		  * Computes a probabilistic upper bound on the error reduction of 
		  * a slice based on an estimate for the mean value in the partition,
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
		upper_bound_error_reduction_mean_estimate_in(
			const unsigned int center, const unsigned int radius,
		  	const std::vector<unsigned int>& partition,
		   	const double average, const double weighted_mean, 
			const unsigned int num_samples_part) const;
		
		/**
		 * Computes a probabilistic upper bound on the error reduction of 
		 * a slice based on the number of vertices inside the partition sampled
		 * in a biased sample. Because this might not biased sampling, just 
		 * returns a very large number.
		 * @param center center
		 * @param radius radius
		 * @param partition partition
		 * @param num_samples_part number of samples used to compute the 
		 * weighted mean.
		 * @throws 
		 * @return upper bound.
		**/
		virtual std::pair<double, double>
		upper_bound_error_reduction_num_samples
			(const unsigned int center, 
			const unsigned int radius,
			const std::vector<unsigned int>& partition,
			const unsigned int num_samples_part,
			const unsigned int total_samples) const
		{
			std::pair<double,double> res;
			res.first = std::numeric_limits<double>::max();
			res.second = 0;
			return res;
		}
		
		virtual double compute_estimate(const double one,
			const double two, const double three) const
		{
			return 0;
		}
};

/**
 * Class that implements the slice tree compression using biased sampling
 * to identify probabilistic good slices.
**/
class SliceTreeBiasSamp: public SliceTreeSamp
{
	public:
		/**
		 * Constructor. Does nothing. 
		 * @param graph graph 
		 * @param max_radius maximum radius for slice tree
		 * @param exhaustive_split consider all splits in all partitions if set
		 * @param delta probability for bounds in error reduction estimates
		 * for slices
		 * @param sampling_rate sampling rate
		 * @return 
		 * @throws 
		**/
		SliceTreeBiasSamp(Graph& _graph, const unsigned int _max_radius, 
			const bool _exhaustive_split, 
			const double _delta, 
			const double _sampling_rate, 
			const double _rho):
			SliceTreeSamp(_graph, _max_radius, _exhaustive_split, 
			_delta, _sampling_rate, _rho)
		{
			graph->set_biased_sampling();
		}
		
		/**
		 * Destructor. Does nothing. 
		 * @param  
		 * @return 
		 * @throws 
		**/
		virtual ~SliceTreeBiasSamp(){;}
	private:
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
		std::pair<double, double>
		upper_bound_error_reduction_num_samples
			(const unsigned int center, const unsigned int radius,
			const std::vector<unsigned int>& partition,
			const unsigned int num_samples_part,
			const unsigned int total_samples) const;
		
		double compute_estimate(const double one,
			const double two, const double three) const;
};

class SliceTreeUnifSamp: public SliceTreeSamp
{
	public:
		/**
		 * @param graph graph 
		 * @param max_radius maximum radius for slice tree
		 * @param exhaustive_split consider all splits in all partitions if set
		 * @param delta probability for bounds in error reduction estimates
		 * for slices
		 * @param sampling_rate sampling_rate
		 * @return 
		 * @throws 
		**/
		SliceTreeUnifSamp(Graph& _graph, const unsigned int _max_radius, 
			const bool _exhaustive_split, 
			const double _delta, const double _sampling_rate,
			const double _rho):
			SliceTreeSamp(_graph, _max_radius, _exhaustive_split, 
			_delta, _sampling_rate,	_rho)
		{
			graph->set_uniform_sampling();
		}

		/**
		 * Destructor. Does nothing. 
		 * @param  
		 * @return 
		 * @throws 
		**/
		virtual ~SliceTreeUnifSamp(){;}
	private:
		/**
		 * Computes a probabilistic upper bound on the error reduction of 
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
			upper_bound_error_reduction_num_samples
			(const unsigned int center, const unsigned int radius,
			const std::vector<unsigned int>& partition,
			const unsigned int num_samples_part,
			const unsigned int total_samples) const;
		
		double compute_estimate(const double one,
			const double two, const double three) const;
};

/**
 * Average linkage tree node
**/
typedef struct ALNode
{
	double average;
	double difference;
	unsigned int size;
	std::vector<unsigned int> partition;

	struct ALNode* left;
	struct ALNode* right;
}al_node_t;

/**
 * Compares two average link nodes
 * By using this function you get an increasing order.
**/
class CompareALNodes
{
	public:
		bool operator()(const al_node_t* n_one, const al_node_t* n_two) const
		{
			return fabs(n_one->difference) > fabs(n_two->difference);
		}
};

/**
 * Class that implements the average linkage compression
**/
class AverageLinkage: public GraphCompressionAlgorithm
{
	public:
		/**
		 * Constructor. Builds the average linkage tree.
		 * @param graph graph 
		 * @param budget budget
		 * @return 
		 * @throws 
		**/
		AverageLinkage(Graph& graph);
		
		/**
		 * Constructor. Does nothing.
		 * @param input_file_name serialized graph with compressed data 
		 * @param graph graph
		 * @return 
		 * @throws 
		**/
		AverageLinkage(const std::string& input_file_name, Graph& graph):
		GraphCompressionAlgorithm(input_file_name, graph){}

		/**
		 * Destructor.
		 * @param
		 * @return 
		 * @throws
		**/
		virtual ~AverageLinkage();
		
		/**
		 * Runs the average linkage compression.
		 * @param
		 * @return 
		 * @throws
		**/
		void compress(const unsigned int budget);

		/**
		 * Decompresses the data from the serialized content of a file
		 * @param input_file_name input file
		 * @return
		 * @throws
		**/
		void decompress(){}

		/**
		 * Prints the average linkage tree
		 * @param st_node parent node
		 * @return
		 * @throws
		**/
		void print();

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
		void write(const std::string& output_file_name)const;
		
		/**
		 * Sets the values as they would be recovered after compression
		 * using average linkage.
		 * @param 
		 * @return
		 * @throws
		**/
		void set_compressed_values();
	private:
		al_node_t* tree;
		double** distance_matrix;
		std::vector< al_node_t* > partitions;
		std::vector<unsigned int> complete_desc_path;
		std::vector<bool> active_partitions;
		unsigned int num_coefficients;
		
		/**
		 * Joins the two last partitions of the complete descending path.
		 * @param
		 * @return 
		 * @throws
		**/
		void join_partitions();

		/**
		 * Tries to perform a basic operation to extend the complete descending path.
		 * @param
		 * @return true in case the extension should continue, false if it is stopped.
		 * @throws
		**/
		bool construct_desc_path();
		
		/**
		 * Computes the difference coefficients for a wavelet
		 * decomposition of the average linkage tree
		 * @param
		 * @return
		 * @throws
		**/
		void compute_difference_coefficients();

		/**
		 * Computes the average coefficients for a wavelet
		 * decomposition of the average linkage tree
		 * @param
		 * @return
		 * @throws
		**/
		void compute_average_coefficients();
		
		/**
		 * Wavelet coefficient pruning of the average linkage tree
		 * @param 
		 * @return
		 * @throws
		**/
		void keep_top_coefficients();
};

/**
 * Wavelet tree node
**/
typedef struct WaveletsNode
{
	double average;
	double difference;
	unsigned int size;
	unsigned int vertex;

	struct WaveletsNode* left;
	struct WaveletsNode* right;
}wavelets_node_t;

/**
 * Compares two wavelets nodes
 * By using this function you get an increasing order.
**/
class CompareWaveletsNodes
{
	public:
		bool operator()(const wavelets_node_t* n_one, 
			const wavelets_node_t* n_two) const
		{
			return fabs(n_one->size * n_one->difference) > 
				fabs(n_two->size * n_two->difference);
		}
};

/**
 * Class that implements the wavelets compression
**/
class Wavelets: public GraphCompressionAlgorithm
{
	public:
		/**
		 * Constructor. Builds the wavelets tree.
		 * @param graph graph 
		 * @param budget budget
		 * @return 
		 * @throws 
		**/
		Wavelets(Graph& graph);
		
		/**
		 * Constructor. Does nothing.
		 * @param input_file_name serialized graph with compressed data 
		 * @param graph graph
		 * @return 
		 * @throws 
		**/
		Wavelets(const std::string& input_file_name, Graph& graph):
		GraphCompressionAlgorithm(input_file_name, graph){}

		/**
		 * Destructor.
		 * @param
		 * @return 
		 * @throws
		**/
		virtual ~Wavelets();
		
		/**
		 * Runs the wavelet compression.
		 * @param
		 * @return 
		 * @throws
		**/
		void compress(const unsigned int budget);

		/**
		 * Decompresses the data from the serialized content of a file
		 * @param input_file_name input file
		 * @return
		 * @throws
		**/
		void decompress(){}

		/**
		 * Prints the wavelet tree
		 * @param
		 * @return
		 * @throws
		**/
		void print();

		/**
		 * Writes the wavelet coefficients of the wavelet compression
		 * to a file.
		 * @param output_file_name ouput file name
		 * @return
		 * @throws
		**/
		void write(const std::string& output_file_name) const;
		
		/**
		 * Sets the values as they would be recovered after compression
		 * using wavelets.
		 * @param 
		 * @return
		 * @throws
		**/
		void set_compressed_values();
	private:
		wavelets_node_t* tree;
		unsigned int num_coefficients;
		
		/**
		 * Computes the average coefficients for a wavelet
		 * tree
		 * @param
		 * @return
		 * @throws
		**/
		void compute_average_coefficients();
		
		/**
		 * Computes the difference coefficients for a wavelet
		 * tree
		 * @param
		 * @return
		 * @throws
		**/
		void compute_difference_coefficients();
		
		/**
		 * Wavelet coefficient pruning
		 * @param 
		 * @return
		 * @throws
		**/
		void keep_top_coefficients();
		
		/**
		 * Builds a wavelet tree from the wavelet nodes
		 * @param wavelets_nodes wavelets nodes
		 * @return
		 * @throws
		**/
		void build_wavelet_tree_recursive(std::vector<wavelets_node_t*>& wavelets_nodes);

		const double compute_sse();

		void build_wavelet_tree();
};

/**
 * Class for graph compression
**/
class GraphCompression
{
	public:
		/**
		 * Compress the graph to the output_file using the available budget.
		 * @param graph_to_compress graph
		 * @param budget budget
		 * @param output_file_name ouput file name
		 * @return a tuple <sse,sse reduction,compresion rate,compression time>
		 * @throws
		 **/
		 static void compress(Graph& graph_to_compress, 
		 	GraphCompressionAlgorithm& algorithm, const unsigned int budget, 
			const std::string& output_file_name);

		/**
		 * Decompress the file to a graph.
		 * @param compressed_file_name compressed file
		 * @return graph
		 * @throws
		 **/
		static void decompress(const std::string& compressed_file_name,
			GraphCompressionAlgorithm& algorithm, Graph& graph);
		
		/**
		 * Returns the sse of the compression
		 * @param 
		 * @return sse
		 * @throws
		 **/
		const static inline double sse()
		{
			return sse_value;
		}

		/**
		 * Returns the sse reduction of the compression
		 * reduction = (sse data - sse compression) / sse data
		 * @param 
		 * @return sse reduction
		 * @throws
		 **/
		const static inline double sse_reduction()
		{
			return sse_reduction_value;
		}

		/**
		 * Returns the compression rate
		 * compression rate = size original data / size compressed data
		 * @param 
		 * @return compression rate
		 * @throws
		 **/
		const static inline double compression_rate()
		{
			return compression_rate_value;
		}

		/**
		 * Returns the budget of the compression
		 * @param 
		 * @return budget
		 * @throws
		 **/
		const static unsigned int budget()
		{
			return budget_value;
		}

		/**
		 * Returns the normalized sse of the compression
		 * normalized sse = sse / sum values squared
		 * @param 
		 * @return normalized sse
		 * @throws
		 **/
		const static double normalized_sse()
		{
			return normalized_sse_value;
		}

		/**
		 * Returns the root mean squared error of the compression
		 * RMSE = sqrt(sse/number of vertices)
		 * @param 
		 * @return root mean squared error
		 * @throws
		 **/
		const static double root_mean_squared_error()
		{
			return root_mean_squared_error_value;
		}

		/**
		 * Returns the maximum pointwise error of the compression
		 * MPE = max difference between value and recovered value in absolute value
		 * @param 
		 * @return maximum pointwise error
		 * @throws
		 **/
		const static double maximum_pointwise_error()
		{
			return maximum_pointwise_error_value;
		}
		
		/**
		 * Returns the peak signal to noise ratio of the compression
		 * PSNR = 20*log_10(max value / RMSE)
		 * @param 
		 * @return maximum pointwise error
		 * @throws
		 **/
		const static double peak_signal_to_noise_ratio()
		{
			return peak_signal_to_noise_ratio_value;
		}
	private:
		static double sse_value;
		static double sse_reduction_value;
		static double compression_rate_value;
		static unsigned int budget_value;
		static double maximum_pointwise_error_value;
		static double peak_signal_to_noise_ratio_value;
		static double root_mean_squared_error_value;
		static double normalized_sse_value;
		static GraphCompressionAlgorithm* compression_algorithm;
		
		/**
		 * Computes several statistics (mostly error measures) for 
		 * the compression
		 * @param graph graph
		 * @return
		 * @throws
		**/
		static void compute_statistics(Graph& graph);
		
		/**
		 * Computes the sse of the compression
		 * @param 
		 * @return sse
		 * @throws
		 **/
		static void compute_sse(Graph& graph);
		
		/**
		 * Computes the sse reduction of the compression
		 * reduction = (sse data - sse compression) / sse data
		 * @param 
		 * @return
		 * @throws
		 **/
		static void compute_sse_reduction(Graph& graph);
		
		/**
		 * Computes the normalized sse of the compression
		 * normalized sse = sse / sum values squared
		 * @param 
		 * @return 
		 * @throws
		 **/
		static void compute_normalized_sse(Graph& graph);

		/**
		 * Computes the compression rate
		 * compression rate = size original data / size compressed data
		 * @param 
		 * @return 
		 * @throws
		 **/
		static void compute_compression_rate(Graph& graph);

		/**
		 * Computes the root mean squared error of the compression
		 * RMSE = sqrt(sse/number of vertices)
		 * @param 
		 * @return 
		 * @throws
		 **/
		static void compute_root_mean_squared_error(Graph& graph);
		
		/**
		 * Computes the maximum pointwise error of the compression
		 * MPE = max difference between value and recovered value in absolute value
		 * @param 
		 * @return 
		 * @throws
		 **/
		static void compute_maximum_pointwise_error(Graph& graph);
		
		/**
		 * Computes the peak signal to noise ratio of the compression
		 * PSNR = 20*log_10(max value / RMSE)
		 * @param 
		 * @return 
		 * @throws
		 **/
		static void compute_peak_signal_to_noise_ratio(Graph& graph);
};

#endif
