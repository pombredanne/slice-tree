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
 *	FILE io.cc: Implementation of I/O operations.
**/

/*my includes*/
#include "io.h"

/*Setting default values for the parameters*/
std::string Parameters::graph_file_name = "";
std::string Parameters::values_file_name = "";
std::string Parameters::compression_algorithm = "";
std::string Parameters::output_file_name = "";
std::string Parameters::partition_sizes_file_name = "";
unsigned int Parameters::budget = 0;
unsigned int Parameters::num_partitions = 0;
double Parameters::sampling_rate = 0;
unsigned int Parameters::num_threads = 1;
double Parameters::delta = 0.9;
unsigned int Parameters::max_radius = USHRT_MAX;
double Parameters::rho = 1;
bool Parameters::directed = false;
bool Parameters::print_tree = false;
bool Parameters::exhaustive_split = false;
std::vector<std::string> Parameters::compression_algorithms;

/**
 * Prints the usage of this program
 * @param
 * @return
 * @throws
**/
void Parameters::print_usage()
{
	std::cout << "Usage: ./graph_compression [OPTIONS]..." << std::endl;
	std::cout << "Compresses the values in a graph using the specified algorithm" << std::endl;
	std::cout << " -g, --graph             input graph file" << std::endl;
	std::cout << " -v, --values            input values file" << std::endl;
	std::cout << " -o, --output            output file name" << std::endl;
	std::cout << " -c, --compression	   compression algorithm" << std::endl;
	std::cout << " -b, --budget            one or more budget values (bytes)" << std::endl;
	std::cout << " -p, --numpart		   one or more numbers of partitions" << std::endl;
	std::cout << " -n, --sampling-rate     sampling rate" << std::endl;
	std::cout << " -s, --partsizes         file with pre-computed partition sizes" << std::endl;
	std::cout << " -t, --numthreads        number of threads to be used" << std::endl;
	std::cout << " -d, --delta             confidence parameter for compression" << std::endl;
	std::cout << " -m, --maxradius         max radius for slice tree compression" << std::endl;
	std::cout << " -r, --rho               approximation constant" << std::endl;
	std::cout << " -e, --directed          if the input graph is directed" << std::endl;
	std::cout << " -i, --print-tree        prints the ST (only when comrression is ST/STUS/STBS) " << std::endl;
	std::cout << " -x, --exhaustive-split  Exhaustively compares slices in all partitions (ST/STBS/STUS only). Significantly slower when set. " << std::endl;
}

/**
 * Reads the input parameters for the compression
 * @param argc size of the command line
 * @param argv string command line
 * @return 
 * @throws invalidParameterSettingException
**/
bool Parameters::read(int argc, char** argv) throw (InvalidParameterSettingException)
{
	InvalidParameterSettingException invalid_parameters;
	
	try
	{
		GetOpt::GetOpt_pp ops(argc, argv);
 
		if (ops >> GetOpt::OptionPresent('e', "directed"))
		{
			directed = true;
		}
		if (ops >> GetOpt::OptionPresent('i', "print-tree"))
		{
			print_tree = true;
		}
		if (ops >> GetOpt::OptionPresent('x', "exhaustive-split"))
		{
			exhaustive_split = true;
		}
		
		if (ops >> GetOpt::OptionPresent('h', "help"))
		{
			print_usage();
			
			return false;
		}

		ops >> GetOpt::Option('g', "graph", graph_file_name, "")
		    >> GetOpt::Option('v', "values", values_file_name, "")
		    >> GetOpt::Option('o', "output", output_file_name, "")
		    >> GetOpt::Option('c', "compression", compression_algorithm, "")
		    >> GetOpt::Option('b', "budget", budget)
		    >> GetOpt::Option('p', "numpart", num_partitions)
		    >> GetOpt::Option('n', "sampling-rate", sampling_rate)
		    >> GetOpt::Option('s', "partsizes", partition_sizes_file_name)
		    >> GetOpt::Option('t', "numthreads", num_threads)
		    >> GetOpt::Option('d', "delta", delta)
		    >> GetOpt::Option('m', "maxradius", max_radius)
		    >> GetOpt::Option('r', "rho", rho);
	}
	catch(GetOpt::GetOptEx& e)
	{
		std::cerr << "Fatal error while parsing the command line parameters!" << std::endl;
		std::cerr << "Try \'./graph_compression --help\' for more information." << std::endl;
		throw invalid_parameters;
	}

	return true;
}

/**
 * Prints the input parameters on the terminal
 * @param
 * @return
 * @throws
**/
void Parameters::print()
{
	std::cout << "graph: " << graph_file_name << std::endl;
	std::cout << "values: " << values_file_name << "\n";
	std::cout << "compression algorithm: " << compression_algorithm << std::endl;
	std::cout << "output: " << output_file_name << std::endl;
	std::cout << "budget: " << budget << std::endl;
	std::cout << "num_partitions: " << num_partitions << std::endl;
}

/**
 * Sets the list of valid compression algorithms, so that 
 * an algorithm given as input can be checked
 * @param algorithms vector with the valid algorithm identifiers
 * @return
 * @throws
 **/
void Parameters::set_compression_algorithms(std::vector<std::string>& algorithms)
{
	compression_algorithms = algorithms;
}

