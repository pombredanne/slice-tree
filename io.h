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
 *	FILE io.h: Definitions of classes related to I/O.
**/

#ifndef IO_H
#define IO_H

/*std includes*/
#include <string>
#include <exception>
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits.h>

/*my includes*/
#include "getopt_pp.h"

/**
 * Parameter setting exception
**/
class InvalidParameterSettingException: public std::exception
{
	virtual const char* what() const throw()
	{
		return "Invalid input parameter setting!";
	}
};

/**
 * Simple handler for the input parameters
**/
class Parameters
{
	public:
		/*Input parameters*/
		static std::string graph_file_name;
		static std::string values_file_name;
		static unsigned int budget;
		static unsigned int num_partitions;
		static std::string compression_algorithm;
		static std::string output_file_name;
		static double sampling_rate;
		static std::string partition_sizes_file_name;
		static unsigned int num_threads;
		static double delta;
		static unsigned int max_radius;
		static double rho;
		static bool directed;
 		static bool print_tree;
		static bool exhaustive_split;
		
		/*List of valid compression algorithms for checking*/
		static std::vector < std::string > compression_algorithms;
		
		/**
		 * Reads the input parameters for the compression
		 * @param argc size of the command line
		 * @param argv string command line
		 * @return 
		 * @throws invalidParameterSettingException
		**/
		static bool read(int argc, char** argv) throw (InvalidParameterSettingException);
		
		/**
		 * Prints the input parameters on the terminal
		 * @param
		 * @return
		 * @throws
		**/
		static void print();

		/**
		 * Sets the list of valid compression algorithms, so that 
		 * an algorithm given as input can be checked
		 * @param algorithms vector with the valid algorithm identifiers
		 * @return
		 * @throws
		 **/
		static void set_compression_algorithms(std::vector<std::string>& algorithms);

		/**
		 * Prints the usage of this program
		 * @param
		 * @return
		 * @throws
		**/
		static void print_usage();
};

#endif
