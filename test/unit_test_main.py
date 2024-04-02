import copy
import numpy as np
import unittest
import sys
sys.path.insert(0, '../')
from graph_table import *


class TestGraphTable(unittest.TestCase):
    def check_orbits(self):
        t_graph = GraphTable(4)
        t_graph.init_single_emitter_graphs(all_connected=True)
        t_graph.generate_orbit_connections()
        self.assertTrue(t_graph.n_orbit == 9)

        t_graph = GraphTable(4)
        t_graph.init_single_emitter_graphs(all_connected=False)
        t_graph.generate_orbit_connections()
        self.assertTrue(t_graph.n_orbit == 13)

if __name__ == '__main__':
    t = TestGraphTable()
    t.check_orbits()

