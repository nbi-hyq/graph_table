import networkx as nx
import numpy as np
import copy
from graph_transformer import measure_single, measure_double_parity
from graph_transformer_with_clifford import measure_double_parity_with_clifford


# locally complement a graph at node m
def local_compl(g, m):
    g_sub = g.subgraph([nb for nb in g.neighbors(m)])
    g_sub_c = nx.complement(g_sub)
    g.remove_edges_from(g_sub.edges())
    g.add_edges_from(g_sub_c.edges())


# convert networkx graph to dense array [num_nodes, nodes, edges], assume less than 256 nodes
# (nx.to_scipy_sparse_array(g): issue with isolated nodes)
def nx_to_array(g):
    num_nodes = g.number_of_nodes()
    arry = np.zeros(2*g.number_of_edges() + num_nodes + 1, dtype=np.uint8)
    arry[0] = num_nodes
    for i_n, node in enumerate(g.nodes):
        arry[i_n+1] = node
    for i, e in enumerate(g.edges):
        arry[2*i+num_nodes+1] = e[0]
        arry[2*i+num_nodes+2] = e[1]
    return arry


# convert dense array representation to networkx graph, assume less than 256 nodes
# (nx.from_scipy_sparse_array(arry): issue with isolated nodes)
def array_to_nx(arry):
    num_nodes = int(arry[0])
    g = nx.Graph()
    g.add_nodes_from([arry[i+1] for i in range(num_nodes)])
    n_edges = (len(arry)-num_nodes-1) // 2
    g.add_edges_from([(arry[2*i+num_nodes+1], arry[2*i+num_nodes+2]) for i in range(n_edges)])
    return g


# print fusion with local Clifford gates that need to be applied afterward to back to a graph state
def print_fusion_with_gates_afterwards(gr_input, fusion_link):
    nx.set_node_attributes(gr_input, '', 'LC')
    measure_double_parity_with_clifford(gr_input, fusion_link[1], fusion_link[2], fusion_link[3])
    l_gates = []
    for n in gr_input.nodes:
        if not gr_input.nodes[n]['LC'] == '':
            l_gates.append(gr_input.nodes[n]['LC'] + '^{-1}_' + str(n))
    print(fusion_link[1:], l_gates)


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

    # find graph in graph dictionary
    def find_graph(self, g):
        h = nx.weisfeiler_lehman_graph_hash(g)
        isomorphic_exists = False
        idx_graph = None
        if h in self.graph_dict.keys():
            for idx_graph in self.graph_dict[h]:
                if nx.is_isomorphic(g, array_to_nx(self.l_graph[idx_graph])):
                    isomorphic_exists = True
                    break
        return h, isomorphic_exists, idx_graph

    # add graph it is not isomorphic to an already existing one
    def add_non_isomorphic(self, g):
        h, isomorphic_exists, idx_graph = self.find_graph(g)
        if isomorphic_exists:
            return isomorphic_exists, idx_graph
        self.l_graph.append(nx_to_array(g))
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

    # initialize from list of given graphs
    def init_from_graphs(self, l_init_graph):
        for g_init in l_init_graph:
            self.add_orbit(g_init)
        self.get_a_graph_orbit()

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
            print(array_to_nx(self.l_graph[idx_gr]).number_of_nodes(), array_to_nx(self.l_graph[idx_gr]).edges)
        print('----------------------------')

    # print the graphs in all orbits
    def print_all_orbits(self, start=0):
        for idx in range(start, self.n_orbit):
            self.print_orbit(idx)

    # create connections between the orbits by measurements
    def generate_orbit_connections(self, measureSingle=False):
        single = ['Z']
        fusion = ['XZZX']
        num_graph_init = self.n_graph
        depth = 0  # depth in the orbit graph
        start = 0  # where the orbits of largest depth start
        while True:
            for i in range(start, num_graph_init):  # only loop over orbits that have been added in previous round
                print(depth, i, self.n_graph)
                gr = array_to_nx(self.l_graph[i])
                _, _, gr_idx = self.find_graph(gr)
                parent_orbit = self.a_graph_orbit[gr_idx]
                if measureSingle:
                    for pauli in single:
                        for n in gr.nodes:
                            g = copy.deepcopy(gr)
                            measure_single(g, n, pauli)
                            self.add_orbit(g, origin=self.l_num_to_orbit[parent_orbit], link=[gr_idx, n, pauli])
                for f in fusion:
                    for n in gr.nodes:
                        for m in gr.nodes:
                            if m < n:
                                g = copy.deepcopy(gr)
                                measure_double_parity(g, n, m, f)
                                self.add_orbit(g, origin=self.l_num_to_orbit[parent_orbit], link=[gr_idx, n, m, f])
            start = num_graph_init
            if num_graph_init == self.n_graph:  # stop if no new graph added
                break
            num_graph_init = self.n_graph
            self.get_a_graph_orbit()
            depth += 1

    # inefficient way of finding connecting path within orbit via BFS (path from gr_idx1 ro gr_idx2)
    def find_path_in_oribit(self, gr_idx1, gr_idx2):
        if gr_idx1 == gr_idx2:
            print('--- no LC ---')
        elif self.a_graph_orbit[gr_idx1] == self.a_graph_orbit[gr_idx2]:
            idx_set = set([gr_idx1])  # to see if graph is there already
            bfs_list = [array_to_nx(self.l_graph[gr_idx1])]
            root_list = [-1]  # graph index in bfs_list on which LC is applied to reach current graph
            lc_list = [-1]  # node number for local complementation
            bfs_pos = 0
            run = True
            # explore with BFS
            while run:
                for n in bfs_list[bfs_pos].nodes:
                    gr_new = copy.deepcopy(bfs_list[bfs_pos])
                    local_compl(gr_new, n)
                    _, _, idx_new = self.find_graph(gr_new)
                    if idx_new not in idx_set:
                        idx_set.add(idx_new)
                        bfs_list.append(copy.deepcopy(gr_new))
                        root_list.append(bfs_pos)
                        lc_list.append(n)
                    if idx_new == gr_idx2:
                        run = False
                        break
                bfs_pos += 1
            # go backwards to find path
            pos = len(root_list) - 1
            path = []
            while pos > 0:
                path.append(lc_list[pos])
                pos = root_list[pos]
            print('LC: ', list(reversed(path)))
        else:
            print('graphs not in same orbit')

    # get the graph index of the graph that is transformed by the measurement
    # link has the form [gr_idx, n_qbit, measurement] or [gr_idx, n_qbit1, n_qbit2, measurement]
    def get_measured_graph_idx(self, link):
        g = array_to_nx(self.l_graph[link[0]])
        if len(link) == 3:
            measure_single(g, link[1], link[2])
        elif len(link) == 4:
            measure_double_parity(g, link[1], link[2], link[3])
        _, _, idx = self.find_graph(g)
        return idx

    # determine way of generating the target graph gr by going backwards in orbit graph
    def back_trace(self, gr):
        _, exist_gr, gr_idx = self.find_graph(gr)
        if exist_gr:
            my_orbit = self.a_graph_orbit[gr_idx]  # orbit index containing the current graph
            print(array_to_nx(self.l_graph[gr_idx]).edges)  # print target graph
            while self.l_num_to_orbit[my_orbit] > 0:
                idx_orbit_entry = self.get_measured_graph_idx(self.l_link_to_orbit[my_orbit])
                self.find_path_in_oribit(idx_orbit_entry, gr_idx)  # LC path from gr_idx to orbit entry point
                print(array_to_nx(self.l_graph[idx_orbit_entry]).edges)  # print graph at orbit entry point
                gr_idx = self.l_link_to_orbit[my_orbit][0]  # move on to graph in upstream orbit from where my_orbit is reached by one fusion
                gr_upstream = array_to_nx(self.l_graph[gr_idx])
                if (self.l_link_to_orbit[my_orbit][3] == 'XZZX'
                        and (self.l_link_to_orbit[my_orbit][1], self.l_link_to_orbit[my_orbit][2]) not in gr_upstream.edges
                        and (self.l_link_to_orbit[my_orbit][2], self.l_link_to_orbit[my_orbit][1]) not in gr_upstream.edges):  # no additional local Clifford gates in this case
                    print(self.l_link_to_orbit[my_orbit][1:])  # print fusion link
                else:
                    print_fusion_with_gates_afterwards(gr_upstream, self.l_link_to_orbit[my_orbit])
                print(array_to_nx(self.l_graph[gr_idx]).edges)  # print graph in upstream orbit
                my_orbit = self.a_graph_orbit[gr_idx]  # move on to upstream orbit containing the new gr_idx graph
            self.find_path_in_oribit(self.l_orbit[my_orbit][0], gr_idx)  # path to branched chain (always at index 0)
            gr_idx_init = self.l_orbit[my_orbit][0]
            print(array_to_nx(self.l_graph[gr_idx_init]).edges)
        else:
            print('graph is not in lookup table')


if __name__ == '__main__':
    import time
    import dill
    load = False
    if load:
        dill.load_module('save_table_14_unconnected_noZ.pkl')
    else:
        t0 = time.time()
        t_graph = GraphTable(14)
        t_graph.init_single_emitter_graphs(all_connected=False)
        t1 = time.time()
        print(t1 - t0)
        t_graph.generate_orbit_connections(measureSingle=False)
        t2 = time.time()
        print(t2 - t1)
        dill.dump_module('save_table_14_unconnected_noZ.pkl')
        t_graph.print_all_orbits(start=t_graph.n_orbit-1)
        print('num_to_orbit frequencies: ', [np.sum(np.array(t_graph.l_num_to_orbit) == i) for i in range(10)])
        print('collisions max: ', t_graph.len_max)

    # make sure the GraphTable object uses the methods defined in this file and not the pickled ones
    from graph_table import *
    t_new = GraphTable(0)
    for v in vars(t_graph).keys():
        exec('t_new.' + v + '= t_graph.' + v)

    # check a few examples
    print('cube: ================================')
    g_cube = nx.Graph()
    g_cube.add_edges_from([(0, 1), (1, 3), (3, 2), (0, 2), (4, 5), (5, 7), (7, 6), (4, 6), (0, 4), (1, 5), (2, 6), (3, 7)])
    t_new.back_trace(g_cube)
    print('stean code: ================================')
    graph = nx.Graph()
    graph.add_edges_from([(0, 1), (1, 3), (3, 2), (0, 2), (4, 5), (4, 6), (0, 4), (1, 5), (2, 6)])
    t_new.back_trace(graph)
    print('5-ring + 2 leave qubits: ================================')
    graph = nx.Graph()
    graph.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (0, 5), (2, 6)])
    t_new.back_trace(graph)
    print('5-ring: ================================')
    graph = nx.Graph()
    graph.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)])
    t_new.back_trace(graph)
    print('6+1 code: ================================')
    graph = nx.Graph()
    graph.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6), (3, 6)])
    t_new.back_trace(graph)
    print('7+1 (FBQC) code from  https://doi.org/10.48550/arXiv.2212.04834 Fig. S3: ================================')
    graph = nx.Graph()
    graph.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6), (1, 6), (3, 6), (4, 6), (2, 7), (5, 7)])
    t_new.back_trace(graph)
    print('8+1 (FBQC) code from  https://doi.org/10.48550/arXiv.2212.04834 Fig. S3: ================================')
    graph = nx.Graph()
    graph.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 0), (0, 7), (3, 7), (4, 7), (2, 8), (5, 8)])
    t_new.back_trace(graph)
    print('9+1 (FBQC) code from  https://doi.org/10.48550/arXiv.2212.04834 Fig. S3: ================================')
    graph = nx.Graph()
    graph.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 0), (0, 7), (3, 7), (4, 7), (2, 8), (5, 8), (1, 9)])
    t_new.back_trace(graph)
    print('10+1 (FBQC) code from  https://doi.org/10.48550/arXiv.2212.04834 Fig. S3: ================================')
    graph = nx.Graph()
    graph.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (7, 9), (0, 8), (0, 9), (1, 8), (6, 9), (8, 10), (9, 10), (2, 10), (5, 10)])
    t_new.back_trace(graph)
    print('cube with one leaf')
    g_cube = nx.Graph()
    g_cube.add_edges_from([(0, 1), (1, 3), (3, 2), (0, 2), (4, 5), (5, 7), (7, 6), (4, 6), (0, 4), (1, 5), (2, 6), (3, 7), (0, 8)])
    t_new.back_trace(g_cube)
    print('from https://doi.org/10.48550/arXiv.2212.04834 Fig. S3: bottom right 7+1 ================================')
    graph = nx.Graph()
    graph.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4), (0, 4), (0, 5), (3, 5), (0, 6), (3, 7)])
    t_new.back_trace(graph)
    print('from https://doi.org/10.48550/arXiv.2212.04834 Fig. S3: bottom right 8+1 ================================')
    graph = nx.Graph()
    graph.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 0), (0, 7), (2, 7), (5, 7), (2, 8), (5, 8)])
    t_new.back_trace(graph)
    print('from https://doi.org/10.48550/arXiv.2212.04834 Fig. S3: bottom right 9+1 ================================')
    graph = nx.Graph()
    graph.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4), (0, 4), (0, 5), (3, 5), (0, 6), (3, 7), (0, 8), (3, 9)])
    t_new.back_trace(graph)

    # repeater graphs
    for n_qubits in range(6, 14, 2):
        print("repeater graph, number of qubits: ", n_qubits)
        half = int(n_qubits / 2)
        graph = nx.Graph()
        graph.add_nodes_from([i for i in range(n_qubits)])
        lEdge = [(i, i + half) for i in range(half)]
        for i in range(half):
            for j in range(i + 1, half):
                lEdge.append((i, j))
        graph.add_edges_from(lEdge)
        t_new.back_trace(graph)

    # crazy graphs
    for n_layer in range(1, 5):
        for size_layer in range(2, 10):
            print("n_layer, size_layer: ", n_layer, size_layer)
            graph = nx.Graph()
            graph.add_nodes_from([i for i in range(n_layer * size_layer + 2)])
            graph.add_edges_from([(i, n_layer * size_layer) for i in range(size_layer)])  # input
            graph.add_edges_from([(i + size_layer * (n_layer - 1), n_layer * size_layer + 1) for i in range(size_layer)])  # output
            for layer in range(n_layer - 1):
                graph.add_edges_from([(i + layer * size_layer, j + (layer + 1) * size_layer) for j in range(size_layer) for i in range(size_layer)])
            t_new.back_trace(graph)
    # wheel graphs
    for size in range(4, 13):
        print("wheel graph size: ", size)
        graph = nx.Graph()
        graph.add_nodes_from([i for i in range(size)])
        graph.add_edges_from([(i, i+1) for i in range(1, size - 1)] + [(1, size - 1)])
        graph.add_edges_from([(0, i) for i in range(1, size)])
        t_new.back_trace(graph)
        nx.draw(graph)

    # the 5 orbits that need 3 fusions
    graph = nx.Graph()
    graph.add_edges_from([(1, 7), (1, 12), (1, 10), (3, 11), (3, 13), (3, 5), (5, 7), (5, 11), (7, 13), (10, 13), (10, 12), (11, 12)])
    t_new.back_trace(graph)
    nx.draw(graph)
    graph = nx.Graph()
    graph.add_edges_from([(1, 4), (1, 6), (1, 7), (1, 9), (1, 10), (3, 9), (3, 4), (3, 7), (4, 11), (6, 7), (6, 11), (9, 10), (10, 11)])
    t_new.back_trace(graph)
    nx.draw(graph)
    graph = nx.Graph()
    graph.add_edges_from([(1, 4), (1, 6), (1, 7), (1, 13), (3, 4), (3, 7), (3, 9), (4, 11), (6, 11), (6, 7), (9, 11), (9, 13), (11, 13)])
    t_new.back_trace(graph)
    nx.draw(graph)
    graph = nx.Graph()
    graph.add_edges_from([(1, 9), (1, 7), (1, 3), (3, 6), (3, 10), (4, 9), (4, 11), (4, 7), (6, 9), (6, 11), (7, 10), (7, 9), (10, 11)])
    t_new.back_trace(graph)
    nx.draw(graph)
    graph = nx.Graph()
    graph.add_edges_from([(1, 9), (1, 3), (1, 7), (3, 10), (3, 6), (4, 9), (4, 11), (4, 7), (6, 9), (6, 11), (7, 10), (10, 11)])
    t_new.back_trace(graph)
    nx.draw(graph)

    # best 1-fusion graph (code)
    graph = nx.Graph()
    graph.add_edges_from([(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 6), (1, 7), (2, 6), (2, 9), (4, 6), (5, 6), (6, 7), (6, 8)])
    t_new.back_trace(graph)
    nx.draw(graph)

    # best 2-fusion graph (code)
    graph = nx.Graph()
    graph.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5), (0, 8), (3, 9), (0, 6), (0, 7), (3, 6), (3, 7)])
    t_new.back_trace(graph)
    nx.draw(graph)

    # count number of connected graph orbits with certain number of nodes
    cnt_connected = [0 for _ in range(14)]
    for i_orbit in range(t_new.n_orbit):
        i_gr = t_new.l_orbit[i_orbit][0]  # take one graph from orbit
        gr0 = array_to_nx(t_new.l_graph[i_gr])
        if len(t_new.l_graph[i_gr]) > 1 and nx.is_connected(gr0):
            cnt_connected[gr0.number_of_nodes() - 1] += 1
    print(cnt_connected)

    # filter graphs that require a certain number of fusions (chose the one with the fewest edges in each orbit)
    num_fusion = 3
    cnt2 = 0
    for i_orbit in range(t_new.n_orbit):
        if t_new.l_num_to_orbit[i_orbit] == num_fusion:
            i_gr = t_new.l_orbit[i_orbit][0]  # pick one graph from orbit for initial filtering
            gr0 = array_to_nx(t_new.l_graph[i_gr])
            if len(t_new.l_graph[i_gr]) > 1 and nx.is_connected(gr0):
                cnt2 += 1
                n_edges_min = np.inf
                for i_gr in t_new.l_orbit[i_orbit]:  # find the graph in the orbit that has the fewest edges
                    gr0 = array_to_nx(t_new.l_graph[i_gr])
                    if gr0.number_of_edges() < n_edges_min:
                        n_edges_min = gr0.number_of_edges()
                        grSelected = copy.deepcopy(gr0)
                print(grSelected.edges)
    print(cnt2)

