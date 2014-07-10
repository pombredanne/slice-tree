#!/usr/bin/python

import sys
import argparse
from graphviz import Digraph

parser = argparse.ArgumentParser(
    description="Generates a GraphViz file from *.graph and *.data files."
    )
parser.add_argument("--data", type=argparse.FileType("r"),
    help="Data input file.")
parser.add_argument("--graph", type=argparse.FileType("r"),
    help="Graph input file.")

if __name__ == "__main__":
    args = parser.parse_args()

    dot = Digraph(
        graph_attr={'size':'64,40'},
        node_attr={'shape':'circle', 'label':'', 'style':'filled',
            'fillcolor':'lightskyblue2', 'nodesep':'1.0', 'ranksep':'1.0'},
        edge_attr={'weight':'1.5'}
        )

    for line in args.data:
        node, value = line.split(",", 1)
        dot.node(node)
    
    for line in args.graph:
        a, b = line.split(",", 1)
        dot.edge(a, b)

    sys.stdout.write(dot.source)
