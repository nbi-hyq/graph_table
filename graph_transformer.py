import networkx as nx
import numpy as np
import random
from functools import reduce

### code from: https://github.com/nbi-hyq/FusionGraphTransformer with local Clifford gates removed

# get node from list (random or deterministic for special neighbor choice)
def get_list_element(l_input, rnd=False):
    if rnd:
        return random.choice(l_input)
    else:
        return max(l_input)


# replace old neighborhood of node n by neighbors in set nb_new
def replace_neigborhood(g, n, nb_new):
    g.remove_edges_from([(n, nb) for nb in g.neighbors(n)])
    for nb in nb_new:
        g.add_edge(n, nb)


# local complementatation (to implement single-qubit measurements or single-qubit gates)
def local_compl(g, node):
    g_sub = g.subgraph([nb for nb in g.neighbors(node)])
    g_sub_c = nx.complement(g_sub)
    g.remove_edges_from(g_sub.edges())
    g.add_edges_from(g_sub_c.edges())


def measure_z(g, node):
    g.remove_node(node)


def measure_y(g, node):
    local_compl(g, node)
    g.remove_node(node)


def measure_y_alternative(g, node):
    nb_node = set([nb for nb in g.neighbors(node)])
    g.remove_node(node)
    nb_n = {}  # dictionary storing original neighborhood of all potentially affected nodes
    for n in nb_node:
        nb_n[n] = set([nb for nb in g.neighbors(n)])
    for n in nb_node:
        replace_neigborhood(g, n, nb_n[n].symmetric_difference(nb_node.difference([n])))


def measure_x(g, node, sp=-1, rnd=False):
    if g.degree[node] == 0:
        g.remove_node(node)
        return
    if sp < 0:
        sp = get_list_element([nb for nb in g.neighbors(node)], rnd=rnd)
    local_compl(g, sp)
    local_compl(g, node)
    g.remove_node(node)
    local_compl(g, sp)


def measure_x_alternative(g, node, sp=-1, rnd=False):
    nb_node = set([nb for nb in g.neighbors(node)])
    g.remove_node(node)
    if not nb_node:
        return
    if sp < 0:
        sp = get_list_element(list(nb_node), rnd=rnd)
    nb_n = {}  # dictionary storing original neighborhood of all potentially affected nodes
    for n in nb_node.union(set([nb for nb in g.neighbors(sp)])):
        nb_n[n] = set([nb for nb in g.neighbors(n)])
    for n in nb_node.difference([sp]):
        if sp not in nb_n[n]:
            replace_neigborhood(g, n, nb_n[n].symmetric_difference(nb_n[sp]).union([sp]))
        else:
            replace_neigborhood(g, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_n[sp], nb_node]).union([sp]))
    for n in nb_n[sp]:
        if n not in nb_node:
            replace_neigborhood(g, n, nb_n[n].symmetric_difference(nb_node))


def transform_xzzx(g_in, qbt_a, qbt_b, a_sp=-1, b_sp=-1):
    nb_a = set([nb for nb in g_in.neighbors(qbt_a)])
    connected = (qbt_b in nb_a)  # connected fusion qubits or not
    if connected:
        nb_a.remove(qbt_b)  # subtract other fusion qubit by default
    nb_b = set([nb for nb in g_in.neighbors(qbt_b)]).difference([qbt_a])  # subtract other fusion qubit by default
    sym_diff_nb_a_nb_b = nb_a.symmetric_difference(nb_b)
    g_in.remove_node(qbt_a)  # remove the fusion qubits after neighborhood is determined
    g_in.remove_node(qbt_b)  # remove the fusion qubits after neighborhood is determined
    nb_n = {}  # dictionary storing original neighborhood of all potentially affected nodes
    if not connected:
        for n in nb_a.union(nb_b):
            nb_n[n] = set([nb for nb in g_in.neighbors(n)])
        for n in nb_a.difference(nb_b):
            replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_b))
        for n in nb_b.difference(nb_a):
            replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_a))
        for n in nb_a.intersection(nb_b):
            replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(sym_diff_nb_a_nb_b))
    else:  # case of connected fusion qubits
        if not (nb_a or nb_b):
            return
        elif (nb_a and not nb_b) or (not nb_a and nb_b):  # lonely qubit
            if not nb_a:  # change meaning if labels A, B if necessary
                nb_h = nb_a
                nb_a = nb_b
                nb_b = nb_h
            if a_sp < 0:
                a_sp = get_list_element(list(nb_a))
            for n in nb_a.union([nb for nb in g_in.neighbors(a_sp)]):
                nb_n[n] = set([nb for nb in g_in.neighbors(n)])
            for n in nb_a.difference([a_sp]):
                if n not in nb_n[a_sp]:
                    replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_n[a_sp]).union([a_sp]))
                else:
                    replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_n[a_sp], nb_a]).union([a_sp]))
            for n in nb_n[a_sp]:
                if n not in nb_a.union(nb_b):
                    replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_a))
        elif not nb_a.difference(nb_b) and not nb_b.difference(nb_a) and nb_a.intersection(nb_b):  # shared neighborhood
            if a_sp < 0:
                a_sp = get_list_element(list(nb_a.intersection(nb_b)))
            for n in nb_a.union([nb for nb in g_in.neighbors(a_sp)]):
                nb_n[n] = set([nb for nb in g_in.neighbors(n)])
            for n in nb_a.intersection(nb_b).difference([a_sp]):
                if n not in nb_n[a_sp]:
                    replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_n[a_sp]).union([a_sp]))
                else:
                    replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_n[a_sp], nb_a]).union([a_sp]))
            for n in nb_n[a_sp]:
                if n not in nb_a.union(nb_b):
                    replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_a))
        elif nb_a.difference(nb_b) and nb_b.difference(nb_a):  # the most general case
            if a_sp < 0:
                a_sp = get_list_element(list(nb_a.difference(nb_b)))
            if b_sp < 0:
                b_sp = get_list_element(list(nb_b.difference(nb_a)))
            for n in nb_a.union(nb_b).union([nb for nb in g_in.neighbors(a_sp)]).union([nb for nb in g_in.neighbors(b_sp)]):
                nb_n[n] = set([nb for nb in g_in.neighbors(n)])
            replace_neigborhood(g_in, a_sp, nb_a.difference([a_sp]))  # other special neighbor done implicitly
            for n in nb_a.difference(nb_b).difference([a_sp]):
                k = nb_n[n].symmetric_difference(nb_n[a_sp]).difference([n, a_sp, b_sp])
                if n not in nb_n[a_sp]:
                    if b_sp not in nb_n[n].symmetric_difference(nb_n[a_sp]):
                        replace_neigborhood(g_in, n, k.union([a_sp]))
                    else:
                        replace_neigborhood(g_in, n, k.symmetric_difference(nb_b.difference([b_sp])).union([a_sp]))
                else:
                    if b_sp not in nb_n[n].symmetric_difference(nb_n[a_sp]):
                        replace_neigborhood(g_in, n, k.symmetric_difference(nb_a.difference([n, a_sp])).union([a_sp]))
                    else:
                        replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [k, nb_b.difference([b_sp]), nb_a.difference([n, a_sp])]).union([a_sp]))
            for n in nb_b.difference(nb_a).difference([b_sp]):
                k = nb_n[n].symmetric_difference(nb_n[b_sp]).difference([n, a_sp, b_sp])
                if n not in nb_n[b_sp]:
                    if a_sp not in nb_n[n].symmetric_difference(nb_n[b_sp]):
                        replace_neigborhood(g_in, n, k.union([b_sp]))
                    else:
                        replace_neigborhood(g_in, n, k.symmetric_difference(nb_a.difference([a_sp])).union([b_sp]))
                else:
                    if a_sp not in nb_n[n].symmetric_difference(nb_n[b_sp]):
                        replace_neigborhood(g_in, n, k.symmetric_difference(nb_b.difference([n, b_sp])).union([b_sp]))
                    else:
                        replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [k, nb_a.difference([a_sp]), nb_b.difference([n, b_sp])]).union([b_sp]))
            for n in nb_a.intersection(nb_b):
                m = reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_n[a_sp], nb_n[b_sp]]).difference([n, a_sp, b_sp])
                if (n not in nb_n[a_sp] and n not in nb_n[b_sp] and b_sp not in nb_n[a_sp]) or (n in nb_n[a_sp] and n in nb_n[b_sp] and b_sp in nb_n[a_sp]):
                    replace_neigborhood(g_in, n, m.union([a_sp, b_sp]))
                elif (n in nb_n[a_sp] and n not in nb_n[b_sp] and b_sp not in nb_n[a_sp]) or (n not in nb_n[a_sp] and n in nb_n[b_sp] and b_sp in nb_n[a_sp]):
                    replace_neigborhood(g_in, n, m.symmetric_difference(nb_a.difference([n, a_sp])).union([a_sp, b_sp]))
                elif (n not in nb_n[a_sp] and n in nb_n[b_sp] and b_sp not in nb_n[a_sp]) or (n in nb_n[a_sp] and n not in nb_n[b_sp] and b_sp in nb_n[a_sp]):
                    replace_neigborhood(g_in, n, m.symmetric_difference(nb_b.difference([n, b_sp])).union([a_sp, b_sp]))
                elif (n not in nb_n[a_sp] and n not in nb_n[b_sp] and b_sp in nb_n[a_sp]) or (n in nb_n[a_sp] and n in nb_n[b_sp] and b_sp not in nb_n[a_sp]):
                    replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [m, nb_a.difference([a_sp]), nb_b.difference([b_sp])]).union([a_sp, b_sp]))
                else:
                    print("error 1")
            for n in nb_n[a_sp].union(nb_n[b_sp]):
                if n not in nb_a.union(nb_b):
                    if n in nb_n[a_sp] and n not in nb_n[b_sp]:
                        replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_a))
                    elif n not in nb_n[a_sp] and n in nb_n[b_sp]:
                        replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_b))
                    else:
                        replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_a, nb_b]))
        elif (nb_a.difference(nb_b) and not nb_b.difference(nb_a)) or (nb_b.difference(nb_a) and not nb_a.difference(nb_b)):  # one neighborhood contained in the other
            if nb_b.difference(nb_a):  # change meaning if labels A, B if necessary
                nb_h = nb_a
                nb_a = nb_b
                nb_b = nb_h
            if a_sp < 0:
                a_sp = get_list_element(list(nb_a.difference(nb_b)))
            if b_sp < 0:
                c_sp = get_list_element(list(nb_b.intersection(nb_a)))
            else:
                c_sp = b_sp  # use specified neighbor (no check if b_sp fulfills conditions)
            for n in nb_a.union(nb_b).union([nb for nb in g_in.neighbors(a_sp)]).union([nb for nb in g_in.neighbors(c_sp)]):
                nb_n[n] = set([nb for nb in g_in.neighbors(n)])
            replace_neigborhood(g_in, c_sp, nb_b.difference([c_sp]))  # other special neighbor done implicitly
            for n in nb_a.difference(nb_b).difference([a_sp]):
                k = nb_n[n].symmetric_difference(nb_n[a_sp]).difference([n, a_sp, c_sp])
                if n not in nb_n[a_sp]:
                    if c_sp not in nb_n[n].symmetric_difference(nb_n[a_sp]):
                        replace_neigborhood(g_in, n, k.union([a_sp]))
                    else:
                        replace_neigborhood(g_in, n, k.symmetric_difference(nb_b.difference([c_sp])).union([a_sp]))
                else:
                    if c_sp not in nb_n[n].symmetric_difference(nb_n[a_sp]):
                        replace_neigborhood(g_in, n, k.symmetric_difference(nb_a.difference(nb_b).difference([n, a_sp])).union([a_sp]))
                    else:
                        replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [k, nb_b.difference([c_sp]), nb_a.difference(nb_b).difference([n, a_sp])]).union([a_sp]))
            for n in nb_a.intersection(nb_b).difference([c_sp]):
                k = nb_n[n].symmetric_difference(nb_n[c_sp]).difference([n, a_sp, c_sp])
                if n not in nb_n[c_sp]:
                    if a_sp not in nb_n[n].symmetric_difference(nb_n[c_sp]):
                        replace_neigborhood(g_in, n, k.union([c_sp]))
                    else:
                        replace_neigborhood(g_in, n, k.symmetric_difference(nb_a.difference(nb_b).difference([a_sp])).union([c_sp]))
                else:
                    if a_sp not in nb_n[n].symmetric_difference(nb_n[c_sp]):
                        replace_neigborhood(g_in, n, k.symmetric_difference(nb_b.difference([n, c_sp])).union([c_sp]))
                    else:
                        replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [k, nb_b.difference([n, c_sp]), nb_a.difference(nb_b).difference([a_sp])]).union([c_sp]))
            for n in nb_n[a_sp].union(nb_n[c_sp]):
                if n not in nb_a.union(nb_b):
                    if n in nb_n[a_sp] and n not in nb_n[c_sp]:
                        replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_a.difference(nb_b)))
                    elif n not in nb_n[a_sp] and n in nb_n[c_sp]:
                        replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_b))
                    else:
                        replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_b, nb_a.difference(nb_b)]))
        else:
            print("error 2")


# XXZZ is identical (independent of whether there is a connection between A, B)
def transform_xxzz(g_in, qbt_a, qbt_b, sp=-1):
    nb_a = set([nb for nb in g_in.neighbors(qbt_a)]).difference([qbt_b])  # subtract other fusion qubit by default
    nb_b = set([nb for nb in g_in.neighbors(qbt_b)]).difference([qbt_a])  # subtract other fusion qubit by default
    sym_diff_nb_a_nb_b = nb_a.symmetric_difference(nb_b)  # The neighbours that only belong to A \Delta B (need for choosing special neighbor)
    g_in.remove_node(qbt_a)  # remove the fusion qubits after neighborhood is determined
    g_in.remove_node(qbt_b)  # remove the fusion qubits after neighborhood is determined
    if not sym_diff_nb_a_nb_b:
        return
    elif sp < 0:
        sp = get_list_element(list(sym_diff_nb_a_nb_b))
    nb_n  = {}  # dictionary storing original neighborhood of all potentially affected nodes
    for n in sym_diff_nb_a_nb_b.union([nb for nb in g_in.neighbors(sp)]):
        nb_n[n] = set([nb for nb in g_in.neighbors(n)])
    for n in sym_diff_nb_a_nb_b.difference([sp]):
        if not n in nb_n[sp]:
            replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_n[sp]).union([sp]))
        else:
            replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_n[sp], sym_diff_nb_a_nb_b]).union([sp]))
    for n in nb_n[sp]:
        if n not in sym_diff_nb_a_nb_b:
            replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_a, nb_b]))


def transform_yzzy(g_in, qbt_a, qbt_b, sp=-1):
    nb_a = set([nb for nb in g_in.neighbors(qbt_a)])
    connected = (qbt_b in nb_a)  # connected fusion qubits or not
    if connected:
        nb_a.remove(qbt_b)  # subtract other fusion qubit by default
    nb_b = set([nb for nb in g_in.neighbors(qbt_b)]).difference([qbt_a])  # subtract other fusion qubit by default
    g_in.remove_node(qbt_a)  # remove the fusion qubits after neighborhood is determined
    g_in.remove_node(qbt_b)  # remove the fusion qubits after neighborhood is determined
    nb_n = {}  # dictionary storing original neighborhood of all potentially affected nodes
    if connected:  # case of connected fusion qubits
        for n in nb_a.union(nb_b):
            nb_n[n] = set([nb for nb in g_in.neighbors(n)])
        for n in nb_a.difference(nb_b):
            replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_a.difference([n])))
        for n in nb_b.difference(nb_a):
            replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_b.difference([n])))
        for n in nb_b.intersection(nb_a):
            replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_b, nb_a]))
    else:  # case of unconnected fusion qubits
        if nb_a.symmetric_difference(nb_b):
            if not nb_b.difference(nb_a):  # swap meaning of A, B
                nb_h = nb_a
                nb_a = nb_b
                nb_b = nb_h
            if sp < 0:
                sp = get_list_element(list(nb_b.difference(nb_a)))
            for n in nb_a.union(nb_b).union([nb for nb in g_in.neighbors(sp)]):
                nb_n[n] = set([nb for nb in g_in.neighbors(n)])
            for n in nb_a.difference(nb_b):
                if not n in nb_n[sp]:
                    replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_n[sp], nb_a.difference([n])]).union([sp]))
                else:
                    replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_b, nb_n[sp].difference([n])]).union([sp]))
            for n in nb_b.difference(nb_a).difference([sp]):
                if not n in nb_n[sp]:
                    replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_n[sp]).union([sp]))
                else:
                    replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_n[sp], nb_a, nb_b]).union([sp]))
            for n in nb_a.intersection(nb_b):
                if not n in nb_n[sp]:
                    replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_a.difference([n])))
                else:
                    replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_b.difference([n])))
            for n in nb_n[sp]:
                if n not in nb_a.union(nb_b):
                    replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_a, nb_b]))
        elif nb_a.intersection(nb_b):
            for n in nb_a.intersection(nb_b):
                nb_n[n] = set([nb for nb in g_in.neighbors(n)])
            for n in nb_a.intersection(nb_b):
                replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_b.difference([n])))


# XYYX is identical (independent of whether there is a connection between A, B)
def transform_xyyx(g_in, qbt_a, qbt_b):
    nb_a = set([nb for nb in g_in.neighbors(qbt_a)]).difference([qbt_b])  # subtract other fusion qubit by default
    nb_b = set([nb for nb in g_in.neighbors(qbt_b)]).difference([qbt_a])  # subtract other fusion qubit by default
    sym_diff_nb_a_nb_b = nb_a.symmetric_difference(nb_b)  # The neighbours that only belong to A \Delta B (need for choosing special neighbor)
    g_in.remove_node(qbt_a)  # remove the fusion qubits after neighborhood is determined
    g_in.remove_node(qbt_b)  # remove the fusion qubits after neighborhood is determined
    nb_n = {}  # dictionary storing original neighborhood of all potentially affected nodes
    for n in sym_diff_nb_a_nb_b:
        nb_n[n] = set([nb for nb in g_in.neighbors(n)])
    for n in sym_diff_nb_a_nb_b:
        replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_a.difference([n]), nb_b.difference([n])]))


def transform_xyyz(g_in, qbt_a, qbt_b, sp=-1):
    nb_a = set([nb for nb in g_in.neighbors(qbt_a)])
    connected = (qbt_b in nb_a)  # connected fusion qubits or not
    if connected:
        nb_a.remove(qbt_b)  # subtract other fusion qubit by default
    nb_b = set([nb for nb in g_in.neighbors(qbt_b)]).difference([qbt_a])  # subtract other fusion qubit by default
    g_in.remove_node(qbt_a)  # remove the fusion qubits after neighborhood is determined
    g_in.remove_node(qbt_b)  # remove the fusion qubits after neighborhood is determined
    nb_n = {}  # dictionary storing original neighborhood of all potentially affected nodes
    if not connected:  # case of unconnected fusion qubits
        for n in nb_a.union(nb_b):
            nb_n[n] = set([nb for nb in g_in.neighbors(n)])
        for n in nb_a.difference(nb_b):
            replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_b))
        for n in nb_b.difference(nb_a):
            replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_b.difference([n]), nb_a]))
        for n in nb_a.intersection(nb_b):
            replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_a.difference([n])))
    else:  # connected fusion qubits
        inter_nb_a_nb_b = nb_a.intersection(nb_b)
        if nb_b.difference(nb_a):
            if sp < 0:
                sp = get_list_element(list(nb_b.difference(nb_a)))
            for n in nb_a.union(nb_b).union([nb for nb in g_in.neighbors(sp)]):
                nb_n[n] = set([nb for nb in g_in.neighbors(n)])
            for n in nb_a.difference(nb_b):
                if not n in nb_n[sp]:
                    replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_a.difference([n])))
                else:
                    replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_b, nb_a.difference([n])]))
            for n in nb_b.difference(nb_a).difference([sp]):
                if not n in nb_n[sp]:
                    replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_n[sp]).union([sp]))
                else:
                    replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_n[sp], nb_b]).union([sp]))
            for n in inter_nb_a_nb_b:
                if not n in nb_n[sp]:
                    replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_n[sp], nb_a.difference([n])]).union([sp]))
                else:
                    replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_a, nb_b, nb_n[sp].difference([n])]).union([sp]))
            for n in nb_n[sp]:
                if n not in nb_a.union(nb_b):
                    replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_b))
        elif inter_nb_a_nb_b:
            if sp < 0:
                sp = get_list_element(list(inter_nb_a_nb_b))
            for n in nb_a.difference(nb_b).union(inter_nb_a_nb_b).union([nb for nb in g_in.neighbors(sp)]):
                nb_n[n] = set([nb for nb in g_in.neighbors(n)])
            for n in nb_a.difference(nb_b):
                if n not in nb_n[sp]:
                    replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_b, nb_a.difference([n])]))
                else:
                    replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_a.difference([n])))
            for n in inter_nb_a_nb_b.difference([sp]):
                if not n in nb_n[sp]:
                    replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_n[sp]).union([sp]))
                else:
                    replace_neigborhood(g_in, n, reduce(lambda a, b: a.symmetric_difference(b), [nb_n[n], nb_n[sp], nb_b]).union([sp]))
            for n in nb_n[sp]:
                if n not in nb_a.union(nb_b):
                    replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_b))
        elif nb_a.difference(nb_b):
            for n in nb_a.difference(nb_b):
                nb_n[n] = set([nb for nb in g_in.neighbors(n)])
            for n in nb_a.difference(nb_b):
                replace_neigborhood(g_in, n, nb_n[n].symmetric_difference(nb_a.difference([n])))


# g: graph, node: measured qubit, pauli: measured Pauli operator (transform if Clifford gate is applied before)
def measure_single(g, node, pauli, method=0, sp=-1, rnd=False):
    pauli_new = pauli
    if pauli_new == 'X':
        if method == 0:
            measure_x(g, node, sp=sp, rnd=rnd)
        else:
            measure_x_alternative(g, node, sp=sp, rnd=rnd)
    elif pauli_new == 'Y':
        if method == 0:
            measure_y(g, node)
        else:
            measure_y_alternative(g, node)
    else:
        measure_z(g, node)


# (fusion success) g: graph, node1, node2: measured qubits, parity: measured double-parity (transform if Clifford gate is applied before)
def measure_double_parity(g, node1, node2, parity):
    parity_new = parity
    if parity_new == 'XZZX':
        transform_xzzx(g, node1, node2)
    elif parity_new == 'XXZZ':
        transform_xxzz(g, node1, node2)
    elif parity_new == 'XYYX':
        transform_xyyx(g, node1, node2)
    elif parity_new == 'YZZY':
        transform_yzzy(g, node1, node2)
    elif parity_new == 'XYYZ':
        transform_xyyz(g, node1, node2)
    elif parity_new == 'YXZY':
        transform_xyyz(g, node2, node1)

