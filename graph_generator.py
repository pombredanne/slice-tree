"""
Copyright (c) 2013, Arlei Silva
All rights reserved.
      
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
	  
Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	        
@author: Arlei Silva (arleilps@gmail.com)
"""

import sys
import getopt
import math
import random
import networkx
from collections import deque

class Graph(object):
    def __init__(self, num_vertices, num_edges, num_partitions, radius, max_average, sse, sse_reduction):
	self.num_vertices = num_vertices
	self.num_edges = num_edges
	self.num_partitions = num_partitions
	self.radius = radius
	self.max_average = max_average
	self.sse_compression = sse - sse_reduction
	self.sse_reduction = sse_reduction
	self.tree = []
    
    def new_node(self, partition):
        node = {}
	node["left"] = -1
	node["right"] = -1
	node["average"] = 0
	node["partition"] = partition
	node["size"] = 0

        return node

    def set_graph(self):
	G =  networkx.barabasi_albert_graph(self.num_vertices, self.num_edges)
        self.edges = []
	self.values = []
	
	for i in range(0, self.num_vertices):
	    self.edges.append([])
	    self.values.append(0)
	
	for e in G.edges(data=True):
	    self.edges[e[0]].append(e[1])
	    self.edges[e[1]].append(e[0])

    def print_tree(self):
        for t in range(0, len(self.tree)):
	    print self.tree[t]

        for t in range(0, len(self.partitions)):
	    print self.partitions[t]

    def set_average_rec(self, node, mu, reduction_slice):
	if(node["partition"] == -1):
#	    print "reduction_slice = %lf\n" % reduction_slice
            size_partition = self.tree[node["left"]]["size"]
	    size_parent = node["size"]
	    size_complement = size_parent - size_partition
	    if size_parent*size_partition > 0:
	        mu_partition = mu + math.sqrt(float(reduction_slice * size_complement) / (size_parent*size_partition))
	    else:
	        mu_partition = 0
            self.tree[node["left"]]["average"] = mu_partition
	    self.set_average_rec(self.tree[node["left"]], self.tree[node["left"]]["average"], reduction_slice)
	    if size_complement > 0:
	        mu_complement = float(mu*size_parent - mu_partition*size_partition) / size_complement
	    else:
	        mu_complement = 0

	    self.tree[node["right"]]["average"] = mu_complement
	    self.set_average_rec(self.tree[node["right"]], self.tree[node["right"]]["average"], reduction_slice)
	else:
	    self.partitions[node["partition"]]["average"] = node["average"]

    def set_averages(self):
	self.tree[0]["average"] = random.uniform(0, self.max_average)
	reduction_slice = float(self.sse_reduction) / (self.num_partitions - 1)
	self.set_average_rec(self.tree[0], self.tree[0]["average"], reduction_slice)
    
    def set_partitions(self):
	self.partition_assignments = []
	self.partitions = []
	self.centers = {}
	self.tree = []
	
	for i in range(0, self.num_vertices):
	    self.partition_assignments.append(0)
	
	partition = {}
	partition["average"] = 0
	partition["center"] = -1
	partition["size"] = self.num_vertices
	partition["node"] = 0
	self.tree.append(self.new_node(0))
	self.tree[0]["size"] = self.num_vertices
	
	self.partitions.append(partition)
	
	for p in range(1, self.num_partitions):
	    partition = self.create_new_partition(p, self.radius)
	    self.partitions.append(partition)
    
    def create_new_partition(self, partition_id, radius):
        partition = {}
	partition["average"] = 0 
	partition["center"] = random.randint(0, self.num_vertices-1)
	partition["node"] = len(self.tree)

	while partition["center"] in self.centers:
	    partition["center"] = random.randint(0, self.num_vertices-1)

	parent = self.partition_assignments[partition["center"]]
	
	self.tree[self.partitions[parent]["node"]]["partition"] = -1
	self.tree[self.partitions[parent]["node"]]["left"] = len(self.tree)
	self.tree[self.partitions[parent]["node"]]["right"] = len(self.tree)+1
	
	self.tree.append(self.new_node(partition_id))
	self.tree.append(self.new_node(parent))

	partition["size"] = self.set_vertices_partition(partition_id, partition["center"], parent, radius)
	self.partitions[parent]["size"] = self.partitions[parent]["size"] - partition["size"]
	self.tree[partition["node"]]["size"] = partition["size"]
	
	self.partitions[parent]["node"] =  self.tree[self.partitions[parent]["node"]]["right"]
	self.tree[self.partitions[parent]["node"]]["size"] = self.partitions[parent]["size"] 

        return partition
    
    def set_vertices_partition(self, partition_id, center, parent, radius):
        distances = {}
        size = 0
	q = deque([center])
	distances[center] = 0
	self.values[center] = 0
	self.partition_assignments[center] = partition_id

	while len(q) > 0:
	    u = q.popleft()
	    
	    for z in self.edges[u]:
	        if z not in distances or distances[z] > distances[u] + 1:
		    if distances[u] + 1 <= radius:
                        distances[z] = distances[u] + 1
			q.append(z)
			if self.partition_assignments[z] == parent:
	                    self.partition_assignments[z] = partition_id
        
	for i in range(0, self.num_vertices):
	    if(self.partition_assignments[i] == partition_id):
	        size = size + 1
	
	return size

    def set_values(self):
	self.values = []
	sse_partition = float(self.sse_compression) / self.num_partitions
	self.set_averages()
#	print "sse_partition = %lf\n" % sse_partition
	for i in range(0, self.num_vertices):
	    self.values.append(0)
        
	for p in range(0, self.num_partitions):
	#    print "partition = %d, average = %lf\n" % (p, self.partitions[p]["average"])
	    size_partition = self.partitions[p]["size"]
	    average = self.partitions[p]["average"]
	    std = math.sqrt(float(sse_partition) / size_partition)
	    for i in range(0, self.num_vertices):
	        if self.partition_assignments[i] == p:
	            self.values[i] = random.gauss(average, std)

    def write_values(self, output_file_name):
        output_file = open(output_file_name, 'w')
        for v in range(0, len(self.edges)):
	    output_file.write(str(v)+","+str(self.values[v])+"\n")
	     
	output_file.close() 
    
    def write_graph(self, output_file_name):
        output_file = open(output_file_name, 'w')
        for v in range(0, len(self.edges)):
	    for u in self.edges[v]:
	        if v > u:
		    output_file.write(str(v)+","+str(u)+"\n")
	     
	output_file.close() 

    def write_statistics(self, output_file_name):
        sse_data = self.compute_sse_data()
	sse = self.compute_sse()
	output_file = open(output_file_name, 'w')
	output_file.write("sse_data = "+str(sse_data)+"\n")
	output_file.write("sse_optimal_partitioning = "+str(sse)+"\n")
#	output_file.write("optimal_sse_reduction = "+str(float(self.sse_data-self.sse)/self.sse_data)+"\n")
	output_file.write("optimal_sse_reduction = "+str(float(sse_data-sse))+"\n")
        output_file.write("partition_id,center,average,actual_average,sse,size\n")
	
	for p in range(0, len(self.partitions)):
	    output_file.write(str(p)+","+str(self.partitions[p]["center"])+","+str(self.partitions[p]["average"])+","+str(self.partitions[p]["actual_average"])+","+str(self.partitions[p]["sse"])+","+str(self.partitions[p]["size"])+"\n")
        
        output_file.write("vertex,partition\n")

	for v in range(0, len(self.partition_assignments)):
	    output_file.write(str(v)+","+str(self.partition_assignments[v])+"\n")
	
	output_file.close() 

    def compute_sse(self):
        sse = 0
	for p in range(0, len(self.partitions)):
	    average = 0
	    n = 0

	    for v in range(0, len(self.partition_assignments)):
	        if self.partition_assignments[v] == p:
		    average = average + self.values[v]
                    n = n + 1
	    
	    if n > 0:
                average = float(average) / n
            
	    sse_partition = 0

	    for v in range(0, len(self.partition_assignments)):
	        if self.partition_assignments[v] == p:
		    sse_partition = sse_partition + math.pow(average - self.values[v], 2)
            
	    self.partitions[p]["sse"] = sse_partition
	    self.partitions[p]["actual_average"] = average
	    sse = sse + sse_partition

        return sse

    def compute_sse_data(self):
        sse = 0
	average = 0

	for v in self.values:
	    average = average + v

	average = float(average) / len(self.values)

#	print "average = %lf" % average

	for v in self.values:
            sse = sse + math.pow(average - v, 2)

	return sse

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg
	    
def main(argv=None):
    if argv is None:
        argv = sys.argv

    #	Parameters:
    #		- output file name	o
    #		- number of vertices v
    #		- number of edges new vertex e (preferential attachment)
    #		- number of partitions	p
    #		- partition radius	r
    #		- max partition average	m
    #		- sse reduction c
    #		- num partitions n
    #

    try:
        try:
            opts, input_files = getopt.getopt(argv[1:], "o:v:e:p:r:m:s:c:h", ["output=","num-vertices=","num-edges=","num-partitions=","radius=","max-average=","sse=","reduction=","help"])
        except getopt.error, msg:
            raise Usage(msg)
 
        output_file_name = ""
	num_vertices = 0
	num_edges = 0
	num_partitions = 0
	radius = 0
	max_average = 0
	sse_reduction = 0
	sse = 0

        for opt,arg in opts:
	    if opt in ('-o', '--output'):
	        output_file_name = arg
	    if opt in ('-v', '--num-vertices'):
	        num_vertices = int(arg)
	    if opt in ('-e', '--num-edges'):
	        num_edges = int(arg)
	    if opt in ('-p', '--num-partitions'):
	        num_partitions = int(arg)
	    if opt in ('-r', '--radius'):
	        radius = int(arg)
	    if opt in ('-m', '--max-average'):
	        max_average = float(arg)
	    if opt in ('-s', '--sse'):
	        sse = float(arg)
	    if opt in ('-c', '--reduction'):
	        sse_reduction = float(arg)
	    if opt in ('-h', '--help'):
	        print "python graph_generator.py [-o <output_file>] [-v <num-vertices>] [-e <num-edges>] [-p <num-partitions] [-r <radius>] [-s <sse>] [-c <reduction>]"
	        sys.exit()
       
        g = Graph(num_vertices, num_edges, num_partitions, radius, max_average, sse, sse_reduction)
	    
	g.set_graph()
	g.write_graph(output_file_name + ".graph")
	g.set_partitions()
	g.set_values()
        g.write_values(output_file_name + ".data")
        g.write_statistics(output_file_name + ".stats")

    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2

if __name__ == "__main__":
    sys.exit(main())
