#!/usr/bin/python

import sys
import argparse

parser = argparse.ArgumentParser(
    description="Generates a GraphViz file from *.graph and *.data files."
    )
parser.add_argument("--data", type=argparse.FileType("r"),
    help="Data input file.")
parser.add_argument("--graph", type=argparse.FileType("r"),
    help="Graph input file.")

gv_header = """digraph {
    graph [size="64,40"]
    node [fillcolor=lightskyblue2 label="" nodesep=1.0 ranksep=1.0 shape=circle style=filled]
    edge [weight=1.5]
"""

gv_footer = "}"

if __name__ == "__main__":
    args = parser.parse_args()

    f = sys.stdout

    try:
        f.write(gv_header)

        # TODO: Escape node names
        for line in args.data:
            node, value = line.strip().split(",", 1)
            f.write("    " + node + "\n")

        for line in args.graph:
            a, b = line.strip().split(",", 1)
            f.write("    " + a + " -> " + b + "\n")

        f.write(gv_footer)
    except OSError:
        pass
