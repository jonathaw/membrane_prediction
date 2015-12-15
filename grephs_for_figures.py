#!/usr/bin/env python3.5
# -*- coding: ascii -*-
# from ProcessEntry import create_topo_entry, process_entry
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib
import networkx as nx
from random import randrange, choice
import operator, sys
matplotlib.rc_context(fname="/home/labs/fleishman/jonathaw/.matplotlib/publishable_matplotlibrc")


class Node:
    def __init__(self, energy, direction, start, end):
        self.energy = energy
        self.direction = direction
        self.start = start
        self.end = end

    def __str__(self):
        return "%-3i ~ %-3i %s %-2i" % (self.start, self.end, self.direction, self.energy)

    def __repr__(self):
        return self.__str__()


def many_wins():
    pass


def draw_graph():
    matplotlib.rcParams['axes.linewidth'] = 1.2
    font0 = FontProperties()
    font = font0.copy()
    font.set_weight('semibold')
    fig = plt.gcf()
    ax = plt.subplot(1, 1, 1)

    G = nx.DiGraph()
    rng = range(60)
    source_node = Node(0, 'src', start=-10, end=-10)
    nodes = [Node(randrange(-10, 5, 1), 'fwd' if choice([True, False]) else 'rev', a, a+20) for a in rng]
    [G.add_node(a, label=str(a)) for a in nodes]
    [G.add_edge(source_node, a, weight=a.energy) for a in nodes]
    [G.add_edge(a, b, weight=b.energy) for a in nodes for b in nodes if b.start >= 1+a.end and a.direction != b.direction]

    pred, dist = nx.bellman_ford(G, source_node)
    sorted_dist = sorted(dist.items(), key=operator.itemgetter(1))
    best_path = [sorted_dist[0][0]]
    while best_path[-1] != source_node:
        best_path.append(pred[best_path[-1]])

    best_path = best_path[::-1]
    best_path_edges = []
    for e in G.edges_iter():
        if e[0] in best_path and e[1] in best_path:
            if e[1] == best_path[best_path.index(e[0])+1]:
                best_path_edges.append(e)
    print(best_path)

    for a in sorted_dist:
        print(a, a[0].direction, best_path[-1].direction)
        if a[0].direction != best_path[-1].direction:
            sec_best_path = [a[0]]
            break
    while sec_best_path[-1] != source_node:
        sec_best_path.append(pred[sec_best_path[-1]])
    sec_best_path = sec_best_path[::-1]
    sec_best_path_edges = []
    for e in G.edges_iter():
        if e[0] in sec_best_path and e[1] in sec_best_path:
            if e[1] == sec_best_path[sec_best_path.index(e[0])+1]:
                sec_best_path_edges.append(e)



    # pos = nx.random_layout(G)
    pos = {n: [n.start, n.energy] for n in nodes}
    pos[source_node] = [0, 0]
    nx.draw_networkx_nodes(G, pos, nodelist=[n for n in nodes if n.direction == 'fwd'], node_color='b', alpha=0.4)
    nx.draw_networkx_nodes(G, pos, nodelist=[n for n in nodes if n.direction == 'rev'], node_color='r', alpha=0.4)
    # nx.draw_networkx_nodes(G, pos, nodelist=best_path, node_color='orange', alpha=0.8)

    nx.draw_networkx_nodes(G, pos, nodelist=[source_node], node_size=400, node_color='g')

    nx.draw_networkx_edges(G, pos, edge_color='gray', alpha=0.2, arrows=True, width=0.8)
    nx.draw_networkx_edges(G, pos, edgelist=best_path_edges, edge_color='black', alpha=0.9, width=2.0)
    nx.draw_networkx_edges(G, pos, edgelist=sec_best_path_edges, edge_color='green', alpha=0.9, width=2.0)

    # nx.draw_networkx_labels(G, pos, c='r', font_size=20, font_family='sans-serif')
    nx.draw_networkx_edge_labels(G, pos, edge_labels={a: a[1].energy for a in best_path_edges+sec_best_path_edges})

    plt.xlim([min(rng)-2, max(rng)+2])
    plt.ylim([min(n.energy for n in nodes)-1, max(n.energy for n in nodes)+1])

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.yticks(fontsize=22)
    plt.xticks(fontsize=22)

    plt.xlabel('Window start position', fontsize=24)
    plt.ylabel(r'$\Delta\Delta G^{insertion}$', fontsize=24)

    fig.set_size_inches(18.5, 10.5)
    plt.savefig('/home/labs/fleishman/jonathaw/graph.png', dpi=100)

    plt.show()


if __name__ == '__main__':
    # many_wins()
    draw_graph()