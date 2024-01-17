import networkx as nx
import numpy as np
import copy
import time
import dill
from graph_transformer import measure_single, measure_double_parity


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
        self.l_link_to_orbit = []  # graph and measurement pattern that links to the orbit
        self.len_max = 0  # keep track of the longest list of non-isomorphic graphs that make hash collisions
        self.per = []  # only used for permuting branches

    # all permutations corresponding to num_elements identical elements filled in boxes (num_separate == num_boxes -1)
    # l_max: should start as num_separate + num_elements
    def generate_permutations(self, num_separate, num_elements, l_max, l_current=[]):
        if num_separate > 0:
            l_current.append(0)
            self.generate_permutations(num_separate - 1, num_elements, l_max, l_current=l_current)
            if l_max == len(l_current):
                self.per.append(tuple(l_current))
            l_current.pop()
        if num_elements > 0:
            l_current.append(1)
            self.generate_permutations(num_separate, num_elements - 1, l_max, l_current=l_current)
            if l_max == len(l_current):
                self.per.append(tuple(l_current))
            l_current.pop()

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
                self.len_max = len(self.graph_dict[h])
        else:
            self.graph_dict[h] = [self.n_graph]
        self.n_graph += 1
        return isomorphic_exists, self.n_graph-1

    # explore the entire orbit starting from g (add all except starting node)
    def explore_orbit_bfs(self, gr):
        l_bfs = [gr]
        bfs_pos = 0
        while bfs_pos < len(l_bfs):
            for n in gr.nodes:
                if [nb for nb in l_bfs[bfs_pos].neighbors(n)]:  # for nodes with no neighbors nothing needs to be done
                    g = copy.deepcopy(l_bfs[bfs_pos])
                    local_compl(g, n)
                    exists_nb, idx_nb = self.add_non_isomorphic(g)
                    if not exists_nb:
                        self.l_orbit[-1].append(idx_nb)
                        l_bfs.append(g)
            bfs_pos += 1

    # add graph if it does not exist and also add the corresponding orbit, origin: number of measurements to reach root
    def add_orbit(self, g, origin=-1, link=[]):
        exists_g, idx_g = self.add_non_isomorphic(g)
        if not exists_g:  # go through all graphs that are locally equivalent
            self.l_orbit.append([idx_g])
            self.n_orbit += 1
            self.l_num_to_orbit.append(origin + 1)  # number of measurements to reach the orbit
            self.l_link_to_orbit.append(link)  # incoming line to orbit [graph, node(s), measurement]
            self.explore_orbit_bfs(g)

    # get the map from graph index to its corresponding orbit index
    def get_a_graph_orbit(self):
        self.a_graph_orbit = np.zeros(self.n_graph, dtype=int)  # array indexed by graph index, giving the orbit
        for i_o, orbit in enumerate(self.l_orbit):
            for g_idx in orbit:
                self.a_graph_orbit[g_idx] = i_o

    # initialize all graphs that a single emitter can make, plus local gates/complementation
    def init_single_emitter_graphs(self, all_connected=False):
        for n_start in range(1, self.num_node_max+1):
            if all_connected:
                i_e_start = 2 ** (n_start - 1) - 1
            else:
                i_e_start = 0
            for i_e in range(i_e_start, 2**(n_start-1)):
                gr = nx.empty_graph(n_start)
                edge_pattern = ('{0:0' + str(n_start - 1) + 'b}').format(i_e)  # connections in the line
                for i_k, k in enumerate(edge_pattern):
                    if k == '1':
                        gr.add_edge(i_k, i_k + 1)
                self.add_orbit(gr)
                # add branches to chain (all connected, unconnected is covered by correspondingly longer chain)
                for n_all in range(n_start, self.num_node_max + 1):
                    gr2 = copy.deepcopy(gr)
                    gr2.add_nodes_from([n for n in range(n_start, n_all)])  # nodes for branch arms
                    self.per = []  # reset
                    self.generate_permutations(n_start - 1, n_all - n_start, n_all - 1, l_current=[])
                    for order in self.per:
                        g = copy.deepcopy(gr2)
                        n_line = 0  # node number in line graph to which stuff is attached
                        n_branch = n_start  # node number of branch node that is attached
                        for k in order:
                            if k == 0:
                                n_line += 1
                            else:
                                g.add_edge(n_line, n_branch)
                                n_branch += 1
                        self.add_orbit(g)
            print(n_start, '#orbits', self.n_orbit, '#graphs', self.n_graph)
        self.get_a_graph_orbit()

    # print the graphs in a certain orbit
    def print_orbit(self, idx_orbit):
        for idx_gr in self.l_orbit[idx_orbit]:
            print(self.l_graph[idx_gr].number_of_nodes(), self.l_graph[idx_gr].edges)
        print('----------------------------')

    # print the graphs in all orbits
    def print_all_orbits(self, start=0):
        for idx in range(start, self.n_orbit):
            self.print_orbit(idx)

    # create connections between the orbits by measurements
    def generate_orbit_connections(self):
        single = ['X', 'Y', 'Z']
        fusion = ['XZZX', 'XXZZ', 'XYYX', 'YZZY', 'XYYZ', 'YXZY']
        num_graph_init = self.n_graph
        depth = 0  # depth in the tablebase
        start = 0  # where the orbits of largest depth start
        while True:
            for i in range(start, num_graph_init):  # only loop over orbits that have been added in previous round
                print(depth, i)
                gr = self.l_graph[i]
                _, gr_idx = self.add_non_isomorphic(gr)  # just for finding graph index of existing graph
                parent_orbit = self.a_graph_orbit[gr_idx]
                for pauli in single:
                    for n in gr.nodes:
                        g = copy.deepcopy(gr)
                        measure_single(g, n, pauli)
                        self.add_orbit(g, origin=self.l_num_to_orbit[parent_orbit], link=[gr_idx, n, pauli])
                for f in fusion:
                    for n in gr.nodes:
                        for m in gr.nodes:
                            if m != n:
                                g = copy.deepcopy(gr)
                                measure_double_parity(g, n, m, f)
                                self.add_orbit(g, origin=self.l_num_to_orbit[parent_orbit], link=[gr_idx, n, m, f])
            start = num_graph_init
            if num_graph_init == self.n_graph:  # stop if no new graph added
                break
            num_graph_init = self.n_graph
            self.get_a_graph_orbit()
            depth += 1

    # determine way of generating g by going backwards in tablebase
    def back_trace(self, gr):
        exist_gr, gr_idx = self.add_non_isomorphic(gr)  # just for finding graph index, TBD: adds graph if non-existing
        if exist_gr:
            parent_orbit = self.a_graph_orbit[gr_idx]
            while self.l_num_to_orbit[parent_orbit] > 0:
                print(self.l_link_to_orbit[parent_orbit])
                gr_idx = self.l_link_to_orbit[parent_orbit][0]
                parent_orbit = self.a_graph_orbit[gr_idx]
            print(self.l_graph[gr_idx].edges)
        else:
            print('graph is not in tablebase')


if __name__ == '__main__':
    load = False
    if load:
        dill.load_module('save_table.pkl')
    else:
        t0 = time.time()
        t_graph = GraphTable(12)
        t_graph.init_single_emitter_graphs(all_connected=False)
        t1 = time.time()
        print(t1 - t0)
        t_graph.generate_orbit_connections()
        t2 = time.time()
        print(t2 - t1)
        dill.dump_module('save_table.pkl')
        t_graph.print_all_orbits(start=t_graph.n_orbit-1)
        print('num_to_orbit frequencies: ', [np.sum(np.array(t_graph.l_num_to_orbit) == i) for i in range(10)])
        print('collisions max: ', t_graph.len_max)

    print('cube: ================================')
    g_cube = nx.Graph()
    g_cube.add_nodes_from([i for i in range(8)])
    g_cube.add_edges_from([(0, 1), (1, 3), (3, 2), (0, 2), (4, 5), (5, 7), (7, 6), (4, 6), (0, 4), (1, 5), (2, 6), (3, 7)])
    t_graph.back_trace(g_cube)
    print('stean code: ================================')
    g_cube = nx.Graph()
    g_cube.add_nodes_from([i for i in range(7)])
    g_cube.add_edges_from([(0, 1), (1, 3), (3, 2), (0, 2), (4, 5), (4, 6), (0, 4), (1, 5), (2, 6)])
    t_graph.back_trace(g_cube)
    print('5-ring + 2 leave qubits: ================================')
    g_cube = nx.Graph()
    g_cube.add_nodes_from([i for i in range(7)])
    g_cube.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (0, 5), (2, 6)])
    t_graph.back_trace(g_cube)
