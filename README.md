# graph_table
This code can be used for computing a lookup table for graph-state generation protocols. It is assumed that a single quantum emitter with a spin can be employed to generate initial graph states with the structure of so-called caterpillar trees. Linear optics type-II fusions are then applied to these initial caterpillar tree graph states to generate a desired target state. After computing th lookup table, a protocol for generating a target graph state with a minimum number of fusions can be found in the lookup table.

# background
- A single quantum emitter one can generate photonic graph states with a caterpillar structure: https://doi.org/10.1103/PhysRevLett.103.113602
- Linear optics type-II fusions and effect on graph states: https://doi.org/10.1103/PhysRevLett.95.010501, https://doi.org/10.48550/arXiv.2405.02414

# source code
- ```graph_table.py``` code to compute the lookup table
- ```graph_table_load.py``` load the lookup table to find some graph-state constructions. To search graph states interactively, load it e.g. in the debug mode.
- ```graph_transformer.py``` graph transformation from https://github.com/nbi-hyq/FusionGraphTransformer with local Clifford gates removed

# example: use the lookup table
Say we want to generate a target graph state with the structure of a cube. Given a list with all its edges, you can find a construction with the following code:

```
graph = nx.Graph()
graph.add_edges_from([(0, 1), (1, 3), (3, 2), (0, 2), (4, 5), (5, 7), (7, 6), (4, 6), (0, 4), (1, 5), (2, 6), (3, 7)])
lookuptable.back_trace(graph)
```

The output will look like:
```
1) gr1 == gr2
2) [27226126, 8, 2, 'XZZX']
3) [(1, 2), (1, 3), (2, 10), (2, 4), (2, 6), (3, 6), (3, 10), (4, 11), (6, 7), (6, 11), (7, 8), (8, 9), (9, 10), (10, 11)]
4) LC-path:  [1, 2, 1]
5) [23535221, 13, 12, 'XZZX']
6) [(1, 2), (1, 4), (1, 12), (1, 6), (2, 3), (3, 4), (4, 11), (6, 7), (6, 11), (7, 8), (8, 9), (9, 10), (10, 13), (11, 12)]
7) gr1 == gr2
8) [4643509, 5, 0, 'XZZX']
9) [(0, 1), (0, 11), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (5, 12), (6, 7), (7, 8), (8, 9), (9, 10), (10, 13)]
10) gr1 == gr2
11) [(0, 1), (0, 11), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (5, 12), (6, 7), (7, 8), (8, 9), (9, 10), (10, 13)]
```

This protocol needs to be read from the bottom to the top. Line 11 is the  graph state that you start with (a caterpillar tree). ```gr1 == gr2``` means that you do not need to apply any local complementations (https://arxiv.org/abs/quant-ph/0602096) and you just keep the graph state as is (line 9). Line 8 then tells you that you need to apply a fusion that measures the parities 'XZ and ZX' (see https://doi.org/10.48550/arXiv.2405.02414) between the qubits 5 and 0 (The number 4643509 tells you where to find the parent graph state in line 9 in the data structure but you can ignore it). ```gr1 == gr2``` in line 7 again says that no local complementations are required before the next fusion, so the graph in line 6 is the graph after the fusion. Line 5 then tells you that you need to apply a fusion between the qubits 12 and 13. Line 4 means that you have to apply some local graph complementations (LC) before the next fusion is performed. The local complementation on the left of the list must be done first, so apply LC(1), LC(2), LC(1). This transforms the graph to the graph shown in line 3. Finally, apply a fusion between qubits 8 and 2 (line 2) and do no more local complementations (```gr1 == gr2``` in line 1). What you get is the cube graph state.





 
