# Slice Tree

An implementation of the Slice Tree algorithm to summarize a network.

# Introduction

In many real-life dynamic graphs, vertex attributes can be represented as
numerical values. For instance, in a traffic network, the average speeds in a
particular location can be modeled as a numerical attribute that changes over
time. Due to the large scale of these attributed graphs, compressing the values
of the attributes is an effective strategy to reduce their processing,
communication, and storage requirements. This research proposes a new graph
compression scheme based on a value-sensitive graph decomposition, which we call
a slice tree. Slices are defined as spherical partitions that capture regions
with homogeneous values in the graph. Values are compressed to the average value
of their respective partitions, minimizing the overall compression error. We
show that our approach outperforms baseline strategies in four real datasets,
achieving error reductions up to 78%.

# Usage

To compile:
```bash
make
```

To generate a synthetic graph with 100000 nodes and 5 edges per node (on
average):
```bash
./graph_generator.py --output synthetic --num-vertices 100000 --num-edges 5 \
    --num-partitions 32 --radius 2 --sse 200000 --reduction 100000
```

To compress the graph with the Slice Tree method using 8 cores:
```bash
./graph_compression --compression ST \
    --values synthetic.data --graph synthetic.graph --partsizes synthetic.index \
    --maxradius 2 --numpart 20 --numthreads 8
```

To compress the graph with the Slice Tree using sampling:
```bash
./graph_compression --compression STBS \
    --values synthetic.data --graph synthetic.graph --partsizes synthetic.index \
    --maxradius 2 --numpart 20 --numthreads 8 --sampling-rate 0.01 --rho 0.9 --delta 0.1
```

To generate a tree, add the `--print-tree` option:
```bash
./graph_compression --compression STBS \
    --values synthetic.data --graph synthetic.graph --partsizes synthetic.index \
    --maxradius 2 --numpart 20 --numthreads 8 --sampling-rate 0.01 --rho 0.9 --delta 0.1 \
    --print-tree > synthetic.tree
```

This can be used to generate an image with [GraphViz][]:

```bash
./viztree synthetic synthetic.data
eog synthetic.svg &
```

[GraphViz]: http://www.graphviz.org
