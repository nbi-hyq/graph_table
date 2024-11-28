# graph_table
This code can be used for using/computing a lookup table for graph-state generation protocols. It is assumed that a single quantum emitter with a spin generates initial graph states with the structure of so-called [caterpillar trees](https://en.wikipedia.org/wiki/Caterpillar_tree). Then, linear optics type-II fusions are applied to these caterpillar tree graph states to generate a desired target state. With the lookup table, you can find a protocol for generating a target graph state with a minimum number of fusions.

# background
- A single quantum emitter can generate photonic graph states with a caterpillar structure: https://doi.org/10.1103/PhysRevLett.103.113602
- Linear optics type-II fusions and effect on graph states: https://doi.org/10.1103/PhysRevLett.95.010501, https://doi.org/10.48550/arXiv.2405.02414

# getting started, source code
- Use ```graph_table_load.py``` to load the lookup table from a ```.pkl```-file and find graph-state constructions. Before running this script, download the ```.pkl```-file that contains the lookup table from https://erda.ku.dk/archives/c34f263db0bbc3649990806548a9833d/published-archive.html and make sure you have the libraries ```networkx, dill, matplotlib, numpy```. To use the lookup table interactively, load it e.g. in the debug mode. Use a system with at least ```32GB``` of RAM.
- ```graph_table.py``` contains code that can be used to compute the lookup table from scratch (takes very long).

# example: use the lookup table
Several examples for lookups of graph state constructions can be found in ```graph_table_load.py```. Let's look at a specific example. Say we want to generate a target graph state with the structure of a cube. Given a list with all its edges and using ```networkx```, you can find a construction with the following code:

```
graph = networkx.Graph()
graph.add_edges_from([(0, 1), (1, 3), (3, 2), (0, 2), (4, 5), (5, 7), (7, 6), (4, 6), (0, 4), (1, 5), (2, 6), (3, 7)])  # add all edges of the cube
t_new.back_trace(graph)  # t_new is how the lookup table object is called in the provided ```.pkl```-file
```

The output will look like:
```
1) [(1, 9), (1, 3), (1, 7), (3, 10), (3, 6), (4, 9), (4, 11), (4, 7), (6, 9), (6, 11), (7, 10), (10, 11)]
2) --- no LC ---
3) [(1, 9), (1, 3), (1, 7), (3, 10), (3, 6), (4, 9), (4, 11), (4, 7), (6, 9), (6, 11), (7, 10), (10, 11)]
4) [8, 2, 'XZZX']
5) [(1, 2), (1, 3), (2, 10), (2, 4), (2, 6), (3, 6), (3, 10), (4, 11), (6, 7), (6, 11), (7, 8), (8, 9), (9, 10), (10, 11)]
6) LC:  [1, 2, 1]
7) [(1, 10), (1, 2), (1, 4), (1, 6), (2, 3), (3, 4), (4, 11), (6, 7), (6, 11), (7, 8), (8, 9), (9, 10), (10, 11)]
8) [13, 12, 'XZZX']
9) [(1, 2), (1, 4), (1, 12), (1, 6), (2, 3), (3, 4), (4, 11), (6, 7), (6, 11), (7, 8), (8, 9), (9, 10), (10, 13), (11, 12)]
10) --- no LC ---
11) [(1, 2), (1, 4), (1, 12), (1, 6), (2, 3), (3, 4), (4, 11), (6, 7), (6, 11), (7, 8), (8, 9), (9, 10), (10, 13), (11, 12)]
12) [5, 0, 'XZZX']
13) [(0, 1), (0, 11), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (5, 12), (6, 7), (7, 8), (8, 9), (9, 10), (10, 13)]
14) --- no LC ---
15) [(0, 1), (0, 11), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (5, 12), (6, 7), (7, 8), (8, 9), (9, 10), (10, 13)]
```

This protocol needs to be read from the bottom to the top. The uneven line numbers represent graph states and the even lines represent fusions or local complementations (see https://arxiv.org/abs/quant-ph/0602096). Line 15 is the graph state with which the physical graph state-generation starts (a caterpillar tree). In l. 14, ```--- no LC ---``` means that you do not need to apply any local complementations and you just keep the graph state as is (so the graphs in l. 13 coincides with the one in l.15). Line 12 says that one needs to apply a fusion (to the state in l.13) that measures the parities ```XZ``` and ```ZX``` (see https://doi.org/10.48550/arXiv.2405.02414) between the qubits 5 and 0. In l.10, ```--- no LC ---``` again says that no local complementations are required before the next fusion, so l.9 equals l.11. Next, one needs to apply a fusion between the qubits 12 and 13 as l.8 indicates. Afterwards, one needs to apply some local graph complementations (LC) before the next fusion is performed (see l.6). The local complementation on the left of the list must be done first, so apply LC(1), LC(2), LC(1). This transforms the graph to the graph shown in l.5. Finally, apply a fusion between qubits 8 and 2 (see l.4) and do no more local complementations (l.2, ```--- no LC ---```). What you get is the cube graph state (see l.1). An illustration is attached below where fusion qubits are highlighted in green.

[cube_construction_all.pdf](https://github.com/user-attachments/files/17916496/cube_construction_all.pdf)

