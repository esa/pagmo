from PyGMO.problem._problem import tsp, tsp_cs, tsp_vrplc, _tsp_encoding
from PyGMO import __extensions__
if __extensions__["gtop"] is True:
    from PyGMO.problem._problem_space import tsp_ds
    tsp_ds.encoding_type = _tsp_encoding

# Renaming and placing the enums
tsp.encoding_type = _tsp_encoding
tsp_vrplc.encoding_type = _tsp_encoding
tsp_cs.encoding_type = _tsp_encoding


def _tsp_ctor(self, weights=[[0, 1, 2], [1, 0, 5], [2, 5, 0]], type="cities"):
    """
        Constructs Travelling Salesman Problem (TSP or ATSP)
        The problem encoding can be of three different types as
        selected by the type kwarg

        1-"cities"
        This encoding represents the ids of the cities visited
        directly in the chromosome. It will
        thus create a constrained problem as only permutation of the
        cities ids are valid (e.g. [0,2,1,5,0] is not
        a valid chromosome)

        2-"randomkeys"
        This encoding, first introduced in the paper
        Bean, J. C. (1994). Genetic algorithms and random keys for
        sequencing and optimization. ORSA journal on computing, 6(2), 154-160.
        It creates a box constrained problem without any constraint.
        It essentially represents the tour as a sequence of doubles bounded
        in [0,1]. The tour is reconstructed by the argsort of the sequence.
        (e.g. [0.34,0.12,0.76,0.03] -> [3,1,0,2])

        3-"full"
        In the full encoding the TSP is represented as a integer linear
        programming problem. The details can be found in
        http://en.wikipedia.org/wiki/Travelling_salesman_problem
        Constructs a Travelling Salesman problem
        (Constrained Integer Single-Objective)

        USAGE: problem.tsp(weights = [0,1,2],[1,0,5],[2,5,0], type="randomkeys")

         * weights: Square matrix with zero diagonal entries containing the cities distances.
         * type: encoding type. One of "cities","randomkeys","full"
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python

    from PyGMO.problem._problem import _tsp_encoding

    def encoding_type(x):
        return {
            "cities": _tsp_encoding.CITIES,
            "randomkeys": _tsp_encoding.RANDOMKEYS,
            "full": _tsp_encoding.FULL
        }[x]

    arg_list = []
    arg_list.append(weights)
    arg_list.append(encoding_type(type))
    self._orig_init(*arg_list)
tsp._orig_init = tsp.__init__
tsp.__init__ = _tsp_ctor


def _tsp_cs_ctor(self, weights=[[0, 1, 2], [1, 0, 5], [2, 5, 0]], values=[1, 1, 1], max_path_length=2, type="cities"):
    """
        Constructs Travelling Salesman Problem City-Selection (TSP-CS)
        The problem encoding can be of three different types as
        selected by the type kwarg

        1-"cities"
        This encoding represents the ids of the cities visited
        directly in the chromosome. It will
        thus create a constrained problem as only permutation of the
        cities ids are valid (e.g. [0,2,1,5,0] is not
        a valid chromosome)

        2-"randomkeys"
        This encoding, first introduced in the paper
        Bean, J. C. (1994). Genetic algorithms and random keys for
        sequencing and optimization. ORSA journal on computing, 6(2), 154-160.
        It creates a box constrained problem without any constraint.
        It essentially represents the tour as a sequence of doubles bounded
        in [0,1]. The tour is reconstructed by the argsort of the sequence.
        (e.g. [0.34,0.12,0.76,0.03] -> [3,1,0,2])

        3-"full"
        In the full encoding the TSP is represented as a integer linear
        programming problem. The details can be found in
        http://en.wikipedia.org/wiki/Travelling_salesman_problem
        Constructs a Travelling Salesman problem
        (Constrained Integer Single-Objective)

        USAGE: problem.tsp_cs(weights=[[0, 1, 2], [1, 0, 5], [2, 5, 0]], values=[1, 1, 1], max_path_length=2, type="cities")

         * weights: Square matrix with zero diagonal entries containing the cities distances.
         * values: The city values.
         * max_path_length: maximum length the salesman can walk
         * type: encoding type. One of "cities","randomkeys","full"
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python

    from PyGMO.problem._problem import _tsp_encoding

    def encoding_type(x):
        return {
            "cities": _tsp_encoding.CITIES,
            "randomkeys": _tsp_encoding.RANDOMKEYS,
            "full": _tsp_encoding.FULL
        }[x]

    arg_list = []
    arg_list.append(weights)
    arg_list.append(values)
    arg_list.append(max_path_length)
    arg_list.append(encoding_type(type))
    self._orig_init(*arg_list)
tsp_cs._orig_init = tsp_cs.__init__
tsp_cs.__init__ = _tsp_cs_ctor


def _tsp_vrplc_ctor(self, weights=[[0, 1, 2], [1, 0, 5], [2, 5, 0]], type="cities", capacity=1.1):
    """
        Constructs Vehicle routing problem with limited capacity.
        This is a variant to the TSP that asks to find n-tours of length
        smaller than the maximum vehicle capacity that visit all cities.
        The objective is to minimize n

        The problem encoding can be of three different types as
        selected by the type kwarg

        1-"cities"
        This encoding represents the ids of the cities visited
        directly in the chromosome. It will
        thus create a constrained problem as only permutation of the
        cities ids are valid (e.g. [0,2,1,5,0] is not
        a valid chromosome)

        2-"randomkeys"
        This encoding, first introduced in the paper
        Bean, J. C. (1994). Genetic algorithms and random keys for
        sequencing and optimization. ORSA journal on computing, 6(2), 154-160.
        It creates a box constrained problem without any constraint.
        It essentially represents the tour as a sequence of doubles bounded
        in [0,1]. The tour is reconstructed by the argsort of the sequence.
        (e.g. [0.34,0.12,0.76,0.03] -> [3,1,0,2])

        3-"full"
        In the full encoding the TSP is represented as a integer linear
        programming problem. The details can be found in
        http://en.wikipedia.org/wiki/Travelling_salesman_problem
        Constructs a Travelling Salesman problem
        (Constrained Integer Single-Objective)

        USAGE: problem.tsp(matrix = [0,1,2],[1,0,5],[2,5,0], type="randomkeys", capacity=1.1)

         * weights: Square matrix with zero diagonal entries containing the cities distances.
         * type: encoding type. One of "cities","randomkeys","full"
         * capacity: maximum vehicle capacity
    """

    from PyGMO.problem._problem import _tsp_encoding

    def encoding_type(x):
        return {
            "cities": _tsp_encoding.CITIES,
            "randomkeys": _tsp_encoding.RANDOMKEYS,
            "full": _tsp_encoding.FULL
        }[x]

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(weights)
    arg_list.append(encoding_type(type))
    arg_list.append(capacity)
    self._orig_init(*arg_list)
tsp_vrplc._orig_init = tsp_vrplc.__init__
tsp_vrplc.__init__ = _tsp_vrplc_ctor


def _plot_tsp(self, x, node_size=10, edge_color='r',
              edge_width=1, bias=None, node_color=None, pos=None):
    """
        Plots a tour represented in the chromosome x
        (using the same encoding of the self object)

        USAGE: problem._plot_tsp(x, node_size=10, edge_color='r',
                edge_width=1, bias=None, node_color=None, pos=None):

         * x:           Crhomosome encoding the city tour.
                        The encoding type used must be the same as that of self
         * node_size:   size of the nodes in the graph visualization
         * edge_color:  size of the edges in the graph visualization
         * edge_width:  width of the edges in the graph visualization
         * bias:        when the graoh node positions are not used,
                        the plot tries to use
                        a spring model to place the nodes. The spring
                        constants depend on this
                        bias parameter
         * node_color:  color of the nodes in the graph visualization
         * pos:         a dictionary containing the node positions
                        (same format as networkx)
    """
    if not (self.verify_x(x) and self.feasibility_x(x)):
        raise Exception("crhomosome is unfeasible")
    from matplotlib import pyplot as plt
    import networkx as nx
    import numpy as np
    from PyGMO.problem import tsp
    fig = plt.gcf()
    axis = plt.gca()

    # We extract few informations on the problem
    weights = self.weights
    n_cities = len(weights[0])
    if self.encoding == _tsp_encoding.RANDOMKEYS:
        edgelist = self.randomkeys2cities(x)
    elif self.encoding == _tsp_encoding.CITIES:
        edgelist = x
    elif self.encoding == _tsp_encoding.FULL:
        edgelist = self.full2cities(x)

    # We construct the list of edges (u,v) containing
    # the indices of the cities visited and we here distinguish between tsp types
    if type(self) == tsp:
        edgelist = [(edgelist[i], edgelist[i + 1]) for i in range(n_cities - 1)] + [(edgelist[-1], edgelist[0])]
    elif type(self) == tsp_cs:
        _, _, id1, id2 = self.find_city_subsequence(x)
        if id1 <= id2:
            edgelist = edgelist[id1:(id2 + 1) % n_cities]
        else:
            edgelist = edgelist[id1:] + edgelist[:id2 + 1]
        edgelist = [(edgelist[i], edgelist[i + 1]) for i in range(len(edgelist) - 1)]
    elif type(self) == tsp_vrplc:
        stl = 0
        chromosome = edgelist
        edgelist = [(chromosome[0], chromosome[1])]
        for i in range(1, n_cities - 1):
            stl += weights[int(chromosome[i])][int(chromosome[i + 1])]
            if stl > self.capacity:
                stl = 0
            else:
                edgelist += [(chromosome[i], chromosome[i + 1])]

    if bias is None:
        bias = max([max(d) for d in weights])

    # We create a networkx graph
    G = nx.Graph()

    # We fill in the vertices
    for i in range(n_cities):
        G.add_node(i)

    # We fill in all the edges
    for i in range(n_cities):
        for j in range(n_cities):
            if i <= j:
                continue
            G.add_edge(i, j, weight=bias / weights[i][j])

    # If cities coordinates are not passed as an input we try to calculate
    # the coordinates for an euclidian TSP (assuming symmetric weights)
    if pos is None:
        # assign the first two nodes: node 0 and node 1, node 0 is
        # chosen to be in the origin
        pos = {0: np.array([0, 0]), 1: np.array([weights[0][1], 0])}
        # algorithm checks during computation of the coordinates if the
        # problem is euclidian
        prob_is_eucl = True
        # we will have to store the first node that is not located in the
        # line constructed by the initial two nodes 0 and 1
        nil_idx = -1
        i = 2
        while (i < n_cities and prob_is_eucl is True):
            # we compute cos(alpha) where alpha is the angle enclosed
            # by the edge (0,1) and (0,i)
            cos_alpha = 0.5 * ((weights[0][i]) ** 2 + (weights[0][1]) ** 2 -
                               (weights[1][i]) ** 2) / (weights[0][i] * weights[0][1])

            if (cos_alpha < -1 or 1 < cos_alpha):
                prob_is_eucl = False
            else:
                # computes one of the two possible positions for node i
                pos[i] = np.array([weights[0][i] * cos_alpha,
                                  weights[0][i] * (1 - cos_alpha ** 2) ** (0.5)])
                omega = 1
                if abs(cos_alpha) != 1:
                    # as soon as one node is not aligned with edge (0,1)
                    # we have to orientate the plot the first node not aligned,
                    # named nil_idx, is chosen to have a positive second
                    # coordinate - every following node is then oriented
                    # accordingly
                    if nil_idx == -1:
                        nil_idx = i
                    elif abs(((pos[i][0] - pos[nil_idx][0]) ** 2 +
                              (pos[i][1] - pos[nil_idx][1]) ** 2) ** (0.5) -
                             weights[i][nil_idx]) > 1e-08 * weights[i][nil_idx]:
                        omega = -1
                pos[i][1] = omega * pos[i][1]  # orient node
                # We have to check the distance to all the previous
                # nodes to decide if the problem is euclidian
                for j in range(2, i):
                    if abs(((pos[i][0] - pos[j][0]) ** 2 +
                           (pos[i][1] - pos[j][1]) ** 2) ** (0.5) -
                       weights[i][j]) > 1e-08 * weights[i][j]:
                        prob_is_eucl = False
            i += 1
        # In case of a non euclidian TSP we create a spring model
        if prob_is_eucl is False:
            pos = nx.layout.spring_layout(G)
    if node_color is None:
        node_color = [0.4] * n_cities

    nx.draw_networkx_nodes(G, pos=pos, node_size=node_size,
                           cmap=plt.cm.Blues, node_color=node_color, ax=axis)
    nx.draw_networkx_edges(G, pos, edgelist=edgelist,
                           width=edge_width, alpha=1, edge_color=edge_color, ax=axis)
    fig.canvas.draw()
    plt.show()
    return pos

tsp.plot = _plot_tsp
tsp_cs.plot = _plot_tsp
tsp_vrplc.plot = _plot_tsp

if __extensions__["gtop"] is True:
    def _tsp_ds_ctor(self, planets, values, max_DV, epochs, type="cities"):
        """
            Constructs Travelling Salesman Problem City-Selection (TSP-CS)
            The problem encoding can be of three different types as
            selected by the type kwarg

            1-"cities"
            This encoding represents the ids of the cities visited
            directly in the chromosome. It will
            thus create a constrained problem as only permutation of the
            cities ids are valid (e.g. [0,2,1,5,0] is not
            a valid chromosome)

            2-"randomkeys"
            This encoding, first introduced in the paper
            Bean, J. C. (1994). Genetic algorithms and random keys for
            sequencing and optimization. ORSA journal on computing, 6(2), 154-160.
            It creates a box constrained problem without any constraint.
            It essentially represents the tour as a sequence of doubles bounded
            in [0,1]. The tour is reconstructed by the argsort of the sequence.
            (e.g. [0.34,0.12,0.76,0.03] -> [3,1,0,2])

            3-"full"
            In the full encoding the TSP is represented as a integer linear
            programming problem. The details can be found in
            http://en.wikipedia.org/wiki/Travelling_salesman_problem
            Constructs a Travelling Salesman problem
            (Constrained Integer Single-Objective)

            USAGE: problem.tsp_cs(planets, values, max_DV, epochs, type="cities"):

             * planets: list of planets
             * values: list of planets values
             * max_DV: maximum DV on-board
             * epochs: list of allowed epochs for the visit (in MJD2000)
             * type: encoding type. One of "cities","randomkeys","full"
        """

        # We construct the arg list for the original constructor exposed by
        # boost_python

        from PyGMO.problem._problem import _tsp_encoding

        def encoding_type(x):
            return {
                "cities": _tsp_encoding.CITIES,
                "randomkeys": _tsp_encoding.RANDOMKEYS,
                "full": _tsp_encoding.FULL
            }[x]

        arg_list = []
        arg_list.append(planets)
        arg_list.append(values)
        arg_list.append(max_DV)
        arg_list.append(epochs)
        arg_list.append(encoding_type(type))
        self._orig_init(*arg_list)
    tsp_ds._orig_init = tsp_ds.__init__
    tsp_ds.__init__ = _tsp_ds_ctor
