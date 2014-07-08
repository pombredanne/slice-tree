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
 *	FILE main.cc: Main method.
**/

/*std includes*/
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

/*my includes*/
#include "io.h"
#include "graph.h"
#include "graph_compression.h"

/**
 * Main method
 * @param argc size of the command line
 * @param argv string containing the command line
 * @return 0 if everything works, 1 otherwise
 * @throws
**/

int main(int argc, char** argv)
{
	
	/*Setting the compression algorithms*/
	std::vector<std::string> compression_algorithms;
	compression_algorithms.push_back("ST");	  //Standard Slice Tree
	compression_algorithms.push_back("STUS"); //Slice Tree with Uniform Sampling
	compression_algorithms.push_back("STBS"); //Slice Tree with Biased Sampling
	compression_algorithms.push_back("STBSR"); //Slice Tree with Biased Sampling and Resampling
	compression_algorithms.push_back("AL");	  //Average Linkage
	compression_algorithms.push_back("WVP");	  //Wavelets with priority vector
	compression_algorithms.push_back("WVB");	  //Wavelets with BFS
	
	Parameters::set_compression_algorithms(compression_algorithms);
		
	double sse;
	double sse_reduction;
	double compression_rate;
	double compression_time;
	double distance_str_time;
	unsigned int num_partitions;
	unsigned int budget;
	double normalized_sse;
	double root_mean_squared_error;
	double maximum_pointwise_error;
	double peak_signal_to_noise_ratio;
	ExecTime* exec_time_compression = new ExecTime();
	GraphCompressionAlgorithm* alg;

	/*Reading the input parameters*/
	if(Parameters::read(argc,argv))
	{
		/*Reading input graph with values*/
		Graph* graph = new Graph(Parameters::graph_file_name, Parameters::values_file_name, Parameters::directed);
		
		if(Parameters::compression_algorithm == "")
		{
			graph->pre_compute_partition_sizes(Parameters::num_threads, 
				Parameters::partition_sizes_file_name,
				Parameters::max_radius);
		}
		else
		{
			if(Parameters::compression_algorithm == "ST" ||
				Parameters::compression_algorithm == "STUS" ||
				Parameters::compression_algorithm == "STBS") 
			{
				graph->read_partition_sizes(Parameters::partition_sizes_file_name,
					Parameters::max_radius);
			}
		}

		/*Performing GraphCompression*/
		if(Parameters::compression_algorithm == "ST")
		{
			/*Standard slice tree*/
			exec_time_compression->start();
			
			if(Parameters::budget > 0)
			{
				alg = new SliceTree(*graph, Parameters::max_radius, Parameters::exhaustive_split); 
				GraphCompression::compress(*graph, *alg, Parameters::budget, 
					Parameters::output_file_name);
			}
			else
			{
				unsigned int budget_from_num_partitions = 
					SliceTree::budget(Parameters::num_partitions, 
					*graph);
				
				alg = new SliceTree(*graph, Parameters::max_radius, Parameters::exhaustive_split); 
				GraphCompression::compress(*graph, *alg, budget_from_num_partitions, 
					Parameters::output_file_name);
			}
			
			exec_time_compression->stop();
			
			if (Parameters::print_tree) {
				SliceTree* st = (SliceTree*) alg; 			
				st->print();
			}
		}
		
		if(Parameters::compression_algorithm == "STUS")
		{
			/*Slice tree with uniform sampling*/
			exec_time_compression->start();
			
			if(Parameters::budget > 0)
			{
				alg = new SliceTreeUnifSamp(*graph, Parameters::max_radius, 
					Parameters::exhaustive_split, 
					Parameters::delta, 
					Parameters::sampling_rate,
					Parameters::rho); 
				GraphCompression::compress(*graph, *alg, Parameters::budget, 
					Parameters::output_file_name);
			}
			else
			{
				unsigned int budget_from_num_partitions = 
					SliceTree::budget(Parameters::num_partitions, 
					*graph);
				
				alg = new SliceTreeUnifSamp(*graph, Parameters::max_radius, Parameters::exhaustive_split,
					Parameters::delta, Parameters::sampling_rate,
					Parameters::rho); 
				GraphCompression::compress(*graph, *alg, budget_from_num_partitions, 
					Parameters::output_file_name);
			}
			
			exec_time_compression->stop();
			if (Parameters::print_tree) {
				SliceTreeUnifSamp* st = (SliceTreeUnifSamp*) alg; 		
				st->print();
			}
		}
		
		if(Parameters::compression_algorithm == "STBS")
		{
			/*Slice tree with biased sampling*/
			exec_time_compression->start();
			
			if(Parameters::budget > 0)
			{
				alg = new SliceTreeBiasSamp(*graph, Parameters::max_radius, Parameters::exhaustive_split,
					Parameters::delta, Parameters::sampling_rate,
					Parameters::rho); 
				GraphCompression::compress(*graph, *alg, Parameters::budget, 
					Parameters::output_file_name);
			}
			else
			{
				unsigned int budget_from_num_partitions = 
					SliceTree::budget(Parameters::num_partitions, 
					*graph);
				alg = new SliceTreeBiasSamp(*graph, Parameters::max_radius, Parameters::exhaustive_split, 
					Parameters::delta, Parameters::sampling_rate, 
					Parameters::rho);
				GraphCompression::compress(*graph, *alg, budget_from_num_partitions, 
					Parameters::output_file_name);
				
			}
			
			exec_time_compression->stop();
			if (Parameters::print_tree) {
				SliceTreeBiasSamp* st = (SliceTreeBiasSamp*) alg; 		
				st->print();
			}
		}
		
		if(Parameters::compression_algorithm == "AL")
		{
			/*Average linkage requires the distance matrix*/
			graph->build_distance_matrix();
			
			exec_time_compression->start();
			
			if(Parameters::budget > 0)
			{
				alg = new AverageLinkage(*graph);
				GraphCompression::compress(*graph, *alg, Parameters::budget, 
					Parameters::output_file_name);
			}
			else
			{
				unsigned int budget_from_num_partitions = 
					SliceTree::budget(Parameters::num_partitions, 
					*graph);
				
				alg = new AverageLinkage(*graph);
				GraphCompression::compress(*graph, *alg, budget_from_num_partitions, 
					Parameters::output_file_name);
			}
			
			exec_time_compression->stop();
		}
		
		if(Parameters::compression_algorithm == "WVP")
		{
			/*Wavelets requires a sorted vector*/
//			graph->build_bfs_vector();
			graph->build_priority_first_vector(0);
			
			exec_time_compression->start();
			
			if(Parameters::budget > 0)
			{
				alg = new Wavelets(*graph);
				GraphCompression::compress(*graph, *alg, Parameters::budget, 
					Parameters::output_file_name);
			}
			else
			{
				unsigned int budget_from_num_partitions = 
					SliceTree::budget(Parameters::num_partitions, 
					*graph);
				alg = new Wavelets(*graph);
				
				GraphCompression::compress(*graph, *alg, budget_from_num_partitions, 
					Parameters::output_file_name);
			}
			
			exec_time_compression->stop();
		}
		
		if(Parameters::compression_algorithm == "WVB")
		{
			/*Wavelets requires a sorted vector*/
			graph->build_bfs_vector();
//			graph->build_priority_first_vector(0);
			
			exec_time_compression->start();
			
			if(Parameters::budget > 0)
			{
				alg = new Wavelets(*graph);
				GraphCompression::compress(*graph, *alg, Parameters::budget, 
					Parameters::output_file_name);
			}
			else
			{
				unsigned int budget_from_num_partitions = 
					SliceTree::budget(Parameters::num_partitions, 
					*graph);
				alg = new Wavelets(*graph);
				
				GraphCompression::compress(*graph, *alg, budget_from_num_partitions, 
					Parameters::output_file_name);
			}
			
			exec_time_compression->stop();
		}

		if(Parameters::compression_algorithm != "")
		{
			sse = GraphCompression::sse();
			sse_reduction = GraphCompression::sse_reduction();
			compression_rate = GraphCompression::compression_rate();
			budget = GraphCompression::budget();
			compression_time = exec_time_compression->get_seconds();
		
			/*Statistics printed as output*/
			std::cout << "budget = " << budget << std::endl;
			std::cout << "sse = " << sse << std::endl;
			std::cout << "sse_reduction = " <<  sse_reduction << std::endl;
			std::cout << "compression_rate = " << compression_rate << std::endl;
			std::cout << "compression_time = " << compression_time << std::endl;
			std::cout << "pruned_slices = " << SliceTreeSamp::pruned() << std::endl;
			
			if(Parameters::compression_algorithm == "STUS" ||
				Parameters::compression_algorithm == "STBS")
			{
				std::cout << "count_bound_1 = " << 
					SliceTreeSamp::count_bound_one() << std::endl;
				
				std::cout << "count_bound_2 = " << 
					SliceTreeSamp::count_bound_two() << std::endl;
				
				std::cout << "count_bound_3 = " << 
					SliceTreeSamp::count_bound_three() << std::endl;
			}
			
			delete alg;
		}

		delete graph;
	}

	delete exec_time_compression;

	return 0;
}


