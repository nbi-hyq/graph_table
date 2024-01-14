import networkx as nx
import numpy as np
import copy
import time
from itertools import permutations


# locally complement a graph at node m
def local_compl(g, m):
    g_sub = g.subgraph([nb for nb in g.neighbors(m)])
    g_sub_c = nx.complement(g_sub)
    g.remove_edges_from(g_sub.edges())
    g.add_edges_from(g_sub_c.edges())


class GraphTable:
    def __init__(self, num_node_max):
        self.num_node_max = num_node_max
        self.graph_dict = {}  # mapping from WL-hash to list of graph indices (in self.l_graph) with that hash
        self.l_adj = []  # measurements linking to other graphs (graph indices), indexed by graph number
        self.l_graph = []  # list of all non-isomorphic graphs, same order as in l_adj
        self.n_graph = 0  # number of all non-isomorphic graphs
        self.l_orbit = []  # list with all the orbits
        self.n_orbit = 0  # number of orbits
        self.a_graph_orbit = None  # array indexed by graph index, giving the orbit
        self.l_num_to_orbit = []  # number of measurements needed to reach orbit, indexed by graph-orbit number
        self.len_max = 0  # keep track of the longest list of non-isomorphic graphs that make hash collisions

    # add graph it is not isomorphic to an already existing one
    def add_non_isomorphic(self, g):
        h = nx.weisfeiler_lehman_graph_hash(g)
        isomorphic_exists = False
        if h in self.graph_dict.keys():
            for idx_graph in self.graph_dict[h]:
                if nx.is_isomorphic(g, self.l_graph[idx_graph]):
                    isomorphic_exists = True
                    break
        if isomorphic_exists:
            return isomorphic_exists, idx_graph
        self.l_graph.append(g)
        if h in self.graph_dict.keys():  # other graph with same hash exists
            self.graph_dict[h].append(self.n_graph)
            if len(self.graph_dict[h]) > self.len_max:
                self.len_max = self.graph_dict[h]
        else:
            self.graph_dict[h] = [self.n_graph]
        self.n_graph += 1
        return isomorphic_exists, self.n_graph-1

    # explore the entire orbit starting from g
    def explore_orbit_recursive(self, gr, idx_g):
        for n in range(gr.number_of_nodes()):
            l_nb = [nb for nb in gr.neighbors(n)]
            if l_nb:  # for nodes with no neighbors nothing needs to be done
                g = copy.deepcopy(gr)  # TBD: needed?
                local_compl(g, n)
                exists_nb, idx_nb = self.add_non_isomorphic(g)
                if not exists_nb and idx_g != idx_nb:  # TBD: do more efficient?
                    self.l_orbit[-1].append(idx_nb)
                    self.explore_orbit_recursive(g, idx_nb)
                #local_compl(g, n)  # restore graph to original

    # initialize all graphs that a single emitter can make, plus local gates/complementation
    def init_single_emitter_graphs(self):
        for n_start in range(1, self.num_node_max+1):
            for i_e in range(2**(n_start-1)):
                for n_all in range(n_start, self.num_node_max + 1):
                    # create graph that can be created by single emitter
                    gr = nx.empty_graph(n_start)
                    gr.add_nodes_from([n for n in range(n_start, n_all)])  # nodes for branch arms
                    edge_pattern = ('{0:0' + str(n_start-1) + 'b}').format(i_e)
                    for i_k, k in enumerate(edge_pattern):
                        if k == '1':
                            gr.add_edge(i_k, i_k+1)
                    # add branches to chain (all connected, unconnected is covered by correspondingly longer chain)
                    l_order = [0 for _ in range(n_start-1)]  # 0 is placeholder between nodes in line
                    l_order.extend([1 for _ in range(n_start, n_all)])
                    per = permutations(l_order)
                    while True:
                        try:
                            order = next(per)
                        except StopIteration:
                            break
                        g = copy.deepcopy(gr)
                        n_line = 0  # node number in line graph to which stuff is attached
                        n_branch = n_start  # node number of branch node that is attached
                        for k in order:
                            if k == 0:
                                n_line += 1
                            else:
                                g.add_edge(n_line, n_branch)
                                n_branch += 1

                        # add the graph if it does not exist
                        exists_g, idx_g = self.add_non_isomorphic(g)

                        # go though all graphs that are locally equivalent
                        if exists_g:
                            continue
                        else:
                            self.l_orbit.append([idx_g])
                            self.n_orbit += 1
                            self.explore_orbit_recursive(g, idx_g)
            print(n_start, '#orbits', self.n_orbit, '#graphs', self.n_graph)

        self.a_graph_orbit = np.zeros(self.n_graph, dtype=int)  # array indexed by graph index, giving the orbit
        for i_o, orbit in enumerate(self.l_orbit):
            for g_idx in orbit:
                self.a_graph_orbit[g_idx] = i_o

    # print the graphs in a certain orbit
    def print_orbit(self, idx_orbit):
        for idx_gr in self.l_orbit[idx_orbit]:
            print(self.l_graph[idx_gr].number_of_nodes(), self.l_graph[idx_gr].edges)
        print('----------------------------')

    # print the graphs in all orbits
    def print_all_orbits(self):
        for idx in range(self.n_orbit):
            self.print_orbit(idx)

    # create connections between the orbits by measurements
    def generate_orbit_connections(self):
        for g in self.l_graph:
            a = None
            # TBD: apply graph tranformations and make links between orbits


if __name__ == '__main__':
    t0 = time.time()
    g_table = GraphTable(8)
    g_table.init_single_emitter_graphs()
    #g_table.print_all_orbits()
    #g_table.print_orbit(10)
    print(time.time()-t0)
    print('coll: ',g_table.len_max)
    g_table.generate_orbit_connections()
