import networkx as nx
import unittest
import sys
sys.path.insert(0, '../')
from graph_table import *


class TestGraphTable(unittest.TestCase):
    def check_orbits(self):
        t_graph = GraphTable(4)
        t_graph.init_single_emitter_graphs(all_connected=True)
        t_graph.generate_orbit_connections(measureSingle=True)
        self.assertTrue(t_graph.n_orbit == 9)

        t_graph = GraphTable(4)
        t_graph.init_single_emitter_graphs(all_connected=True)
        t_graph.generate_orbit_connections(measureSingle=False)
        self.assertTrue(t_graph.n_orbit == 7)

        t_graph = GraphTable(4)
        t_graph.init_single_emitter_graphs(all_connected=False)
        t_graph.generate_orbit_connections(measureSingle=True)
        self.assertTrue(t_graph.n_orbit == 13)

        t_graph = GraphTable(4)
        t_graph.init_single_emitter_graphs(all_connected=False)
        t_graph.generate_orbit_connections(measureSingle=False)
        self.assertTrue(t_graph.n_orbit == 13)

    def test_graph_compression(self):
        for n in [4, 8, 11, 17]:
            for p in [0.1, 0.3, 0.5, 0.7, 0.9]:
                for _ in range(10):
                    g = nx.erdos_renyi_graph(n, p)
                    a = nx_to_array(g)
                    g_recover = array_to_nx(a)
                    self.assertTrue(nx.is_isomorphic(g, g_recover))

    # test assumption made about nx-implementation
    def test_nx(self):
        g1 = nx.Graph()
        g1.add_nodes_from([0, 1, 2])
        g1.add_edges_from([(0, 1), (1, 2)])
        g2 = nx.Graph()
        g2.add_nodes_from([1, 2, 3])
        g2.add_edges_from([(1, 2), (3, 2)])
        self.assertTrue(nx.is_isomorphic(g1, g2))
        self.assertTrue(nx.weisfeiler_lehman_graph_hash(g1) == nx.weisfeiler_lehman_graph_hash(g2))

        g1 = nx.Graph()
        g1.add_nodes_from([0, 1, 2, 3])
        g2 = nx.Graph()
        g2.add_nodes_from([1, 2, 3])
        self.assertFalse(nx.is_isomorphic(g1, g2))
        self.assertFalse(nx.weisfeiler_lehman_graph_hash(g1) == nx.weisfeiler_lehman_graph_hash(g2))


if __name__ == '__main__':
    t = TestGraphTable()
    t.check_orbits()
    t.test_graph_compression()
    t.test_nx()

