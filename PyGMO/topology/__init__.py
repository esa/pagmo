# -*- coding: utf-8 -*-
from PyGMO.topology._topology import *


# Some methods added to interface to networkx
def _to_networkx(self):
    """
    Export topology as a networkx DiGraph.
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError('Could not import the networkx module.')
    retval = nx.DiGraph()
    for i in range(self.number_of_vertices):
        if self.get_num_adjacent_vertices(i):
            retval.add_edges_from([(i, n)
                                  for n in self.get_adjacent_vertices(i)])
        else:
            retval.add_node(i)
    return retval
_topology._base.to_networkx = _to_networkx


def _draw(
        self,
        layout='spring',
        n_color='blue',
        n_size=15,
        n_alpha=0.5,
        e_alpha=0.1,
        e_arrows=False,
        scale_by_degree=False):
    """
    Draw topology using the draw() command from networkx.

    USAGE: topology.draw(layout = "spring", n_size = 15, scale_by_degree = False, n_color = 'blue', n_alpha = 0.5, e_alpha = 0.1, e_arrows=False)

    * layout: Network layout. Can be 'spring' or 'circular'.
    * n_size: The size of nodes. Becomes scaling factor when scale_by_degree=True.
    * scale_by_degree: When True, nodes will be sized proportional to their degree.
    * n_color: Node color.
    * n_alpha: Transparency of nodes. Takes value between 0 and 1.
    * e_elpha: Transparency of edges. Takes value between 0 and 1.
    * e_arrows: Plots arrows on the edges for directed graphs
    * scale_by_degree: When True, nodes will be sized proportional to their degree.
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError('Could not import the networkx module.')
    try:
        import matplotlib.pyplot as pl
    except ImportError:
        raise ImportError('Could not improt the MatPlotLib module.')
    if self.number_of_vertices == 0 or self.number_of_vertices == 1:
        raise ValueError(
            'Cannot draw topology with one single vertex or less.')

    G = self.to_networkx()
    node_sizes = list(range(nx.number_of_nodes(G)))
    for i in range(nx.number_of_nodes(G)):
        if scale_by_degree:
            node_sizes[i] = nx.degree(G, i) * n_size
        else:
            node_sizes[i] = n_size

    if layout == "spring":
        pos = nx.spring_layout(self.to_networkx())
    if layout == "circular":
        pos = nx.circular_layout(self.to_networkx())

    pl.figure()
    nx.draw_networkx_edges(
        self.to_networkx(), pos, alpha=e_alpha, arrows=e_arrows)
    nx.draw_networkx_nodes(
        self.to_networkx(),
        pos,
        node_size=node_sizes,
        node_color=n_color,
        alpha=n_alpha)
    pl.axis('off')
    pl.show()
_topology._base.draw = _draw


def _draw_degree_distribution(self, style='.'):
    """
    Plot the degree dstribution on a loglog scale

    USAGE topology.draw_degree_distribution(style = 'r.')

    style: MatPlotLib line style.
    """
    try:
        import matplotlib.pyplot as pl
    except ImportError:
        raise ImportError('Could not improt the MatPlotLib module.')
    dd = self.get_degree_distribution()
    pl.figure()
    pl.loglog(dd, style)
    pl.xlabel('k - number of links')
    pl.ylabel('p(k) - probability of k')
    pl.show()
_topology._base.draw_degree_distribution = _draw_degree_distribution

# Re-exposing the constructors


def _watts_strogatz_ctor(self, k=10, p=0.1, nodes=0):
    """
    Constructs Watts-Strogatz Topology:

    USAGE: topology.watts-strogatz(k = 10, p = 0.1, nodes=0)

    * k: number of neighbours
    * p: rewiring probability [0,1]
    * nodes: number of nodes
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(k)
    arg_list.append(p)
    arg_list.append(nodes)
    self._orig_init(*arg_list)
watts_strogatz._orig_init = watts_strogatz.__init__
watts_strogatz.__init__ = _watts_strogatz_ctor


def _erdos_renyi_ctor(self, p=0.1, nodes=0):
    """
    Constructs an Erdos-Renyi (random) topology:

    USAGE: topology.erdos_renyi(p = 0.1, nodes=0)

    * p: wiring probability [0,1]
    * nodes: number of nodes
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(p)
    self._orig_init(*arg_list)
    for i in range(nodes):
        self.push_back()
erdos_renyi._orig_init = erdos_renyi.__init__
erdos_renyi.__init__ = _erdos_renyi_ctor


def _barabasi_albert_ctor(self, m0=3, m=3, nodes=0):
    """
    Constructs an Barabasi-Albert topology:

    USAGE: topology.barabasi_albert(m0 = 3, m = 3, nodes = 0)

    * m0: kernel size
    * m: number of random connections to be established when a new node is added.
    * nodes: number of nodes
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(m0)
    arg_list.append(m)
    self._orig_init(*arg_list)
    for i in range(nodes):
        self.push_back()
barabasi_albert._orig_init = barabasi_albert.__init__
barabasi_albert.__init__ = _barabasi_albert_ctor


def _clustered_ba_ctor(self, m0=3, m=3, p=0.5, nodes=0):
    """
    Constructs a Clustered Barabasi-Albert topology:

    USAGE: topology.clustered_ba(mm0 = 3, m=3, p=0.5, nodes = 0)

    * m0: kernel size
    * m: nnumber of random connections to be established when a new node is added
    * p: probability that a connection is established between two nodes that are adjacent to a new node
    * nodes: number of nodes
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(m0)
    arg_list.append(m)
    arg_list.append(p)
    self._orig_init(*arg_list)
    for i in range(nodes):
        self.push_back()
clustered_ba._orig_init = clustered_ba.__init__
clustered_ba.__init__ = _clustered_ba_ctor


def _ageing_clustered_ba_ctor(self, m0=3, m=3, p=0.5, a=1000, nodes=0):
    """
    Constructs a Clustered Barabási-Albert with Ageing vertices graph topology.:

    USAGE: topology.clustered_ba(m0 = 3, m=3, p=0.5, a=1000, nodes = 0)

    * m0: kernel size
    * m: number of random connections to be established when a new node is added
    * p: probability that a connection is established between two nodes that are adjacent to a new node
    * a: 'age' at which a node ceases to make new connections.
    * nodes: number of nodes
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(m0)
    arg_list.append(m)
    arg_list.append(p)
    arg_list.append(a)
    self._orig_init(*arg_list)
    for i in range(nodes):
        self.push_back()
ageing_clustered_ba._orig_init = ageing_clustered_ba.__init__
ageing_clustered_ba.__init__ = _ageing_clustered_ba_ctor


def _fully_connected_ctor(self, nodes=0):
    """
    Constructs a fully connected topology:

    USAGE: topology.fully_connected(nodes = 0)

    * nodes: number of nodes
    """
    # We set the defaults or the kwargs
    arg_list = []
    self._orig_init(*arg_list)
    for i in range(nodes):
        self.push_back()
fully_connected._orig_init = fully_connected.__init__
fully_connected.__init__ = _fully_connected_ctor


def _hypercube_ctor(self, nodes=0):
    """
    Constructs a hypercube topology:

    USAGE: topology.hypercube(nodes = 0)

    * nodes: number of nodes
    """
    # We set the defaults or the kwargs
    arg_list = []
    self._orig_init(*arg_list)
    for i in range(nodes):
        self.push_back()
hypercube._orig_init = hypercube.__init__
hypercube.__init__ = _hypercube_ctor


def _one_way_ring_ctor(self, nodes=0):
    """
    Constructs a one_way_ring topology:

    USAGE: topology.one_way_ring(nodes = 0)

    * nodes: number of nodes
    """
    # We set the defaults or the kwargs
    arg_list = []
    self._orig_init(*arg_list)
    for i in range(nodes):
        self.push_back()
one_way_ring._orig_init = one_way_ring.__init__
one_way_ring.__init__ = _one_way_ring_ctor


def _pan_ctor(self, nodes=0):
    """
    Constructs a pan topology: this is, essentially, a ring with a sink node. This topology was 'invented' to support a local
    optimization algorithm in the sink node

    USAGE: topology.pan(nodes = 0)

    * nodes: number of nodes
    """
    # We set the defaults or the kwargs
    arg_list = []
    self._orig_init(*arg_list)
    for i in range(nodes):
        self.push_back()
pan._orig_init = pan.__init__
pan.__init__ = _pan_ctor


def _rim_ctor(self, nodes=0):
    """
    Constructs a rim topology: this is, essentially, a ring + one node that connects all ring nodes bi-directionally

    USAGE: topology.rim(nodes = 0)

    * nodes: number of nodes
    """
    # We set the defaults or the kwargs
    arg_list = []
    self._orig_init(*arg_list)
    for i in range(nodes):
        self.push_back()
rim._orig_init = rim.__init__
rim.__init__ = _rim_ctor


def _ring_ctor(self, nodes=0):
    """
    Constructs a ring topology

    USAGE: topology.ring(nodes = 0)

    * nodes: number of nodes
    """
    # We set the defaults or the kwargs
    arg_list = []
    self._orig_init(*arg_list)
    for i in range(nodes):
        self.push_back()
ring._orig_init = ring.__init__
ring.__init__ = _ring_ctor


def _unconnected_ctor(self, nodes=0):
    """
    Constructs an unconnected topology

    USAGE: topology.unconnected(nodes = 0)

    * nodes: number of nodes
    """
    # We set the defaults or the kwargs
    arg_list = []
    self._orig_init(*arg_list)
    for i in range(nodes):
        self.push_back()
unconnected._orig_init = unconnected.__init__
unconnected.__init__ = _unconnected_ctor


def _custom_add_edge(self, n, m, weight=1.0):
    """
    Adds and edge to a custom topology with an optional weight argument.

    USAGE: top.add_edge(n, m, weight=1.0)
    * weight: Edge weight between vertices v2 and v2
    """

    if not isinstance(n, int) or not isinstance(m, int):
        raise TypeError("Indices n and m must be integers.")
    if not isinstance(weight, float):
        raise TypeError("Migration probability must be a float.")

    arg_list = []
    arg_list.append(n)
    arg_list.append(m)
    arg_list.append(weight)
    self._orig_add_edge(*arg_list)
custom._orig_add_edge = custom.add_edge
custom.add_edge = _custom_add_edge
