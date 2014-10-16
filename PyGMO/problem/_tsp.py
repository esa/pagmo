from PyGMO.problem._problem import tsp,_tsp_encoding


# Renaming and placing the enums
tsp.encoding = _tsp_encoding

def _tsp_ctor(self, weights=[[0, 1, 2], [1, 0, 5], [2, 5, 0]], type="full"):
    """
		Constructs Travelling Salesman Problem (TSP or ATSP)
		The problem encoding can be of three different types as selcted by the type kwarg
	
		1-"cities"
		This encoding represents the ids of the cities visited directly in the chromosome. It will
		thus create a constrained problem as only permutation of the cities ids are valid (e.g. [0,2,1,5,0] is not
		a valid chromosome)
	
		2-"randomkeys"
		This encoding, first introduced in the paper
		Bean, J. C. (1994). Genetic algorithms and random keys for sequencing and optimization. ORSA journal on computing, 6(2), 154-160.
		creates a box constrained problem without any constraint. It essentially represents the tour as a sequence of doubles bounded in [0,1].
		The tour is reconstructed by the argsort of the sequence. (e.g. [0.34,0.12,0.76,0.03] -> [3,1,0,2])
	
		3-"full"
		In the full encoding the TSP is represented as a integer linear programming problem. The details can be found in
		http://en.wikipedia.org/wiki/Travelling_salesman_problem#Integer_linear_programming_formulation
		Constructs a Travelling Salesman problem (Constrained Integer Single-Objective)

		USAGE: problem.tsp(matrix = [0,1,2],[1,0,5],[2,5,0], type="randomkeys")

		 * weights: Square matrix with zero diagonal entries containing the cities distances.
		 * type: encoding type. One of "cities","randomkeys","full"
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(weights)
    if type == "full":
    	encoding_type = tsp.encoding.FULL
    elif type== "randomkeys":
    	encoding_type = tsp.encoding.RANDOMKEYS
    elif type== "cities":
    	encoding_type = tsp.encoding.CITIES
    else:
    	raise ValueError("Unrecognized encoding type")

    arg_list.append(encoding_type)
    self._orig_init(*arg_list)
tsp._orig_init = tsp.__init__
tsp.__init__ = _tsp_ctor

def _plot_tsp(self,x,pos=None,node_size=10,edge_color='r',edge_width=1,bias=None):
	if not (self.verify_x(x) and self.feasibility_x(x)):
		raise Exception("crhomosome is unfeasible")
	from matplotlib import pyplot as plt
	import networkx as nx

	fig = plt.gcf()

	# We extract few informations on the problem
	weights = self.weights
	n_cities = len(weights[0])
	edgelist = self.randomkeys2cities(x)
	edgelist = [(edgelist[i],edgelist[i+1]) for i in range(n_cities-1)] + [(edgelist[-1],edgelist[0])]
	if bias==None:
		bias = max([max(d) for d in weights])

	# We create a networkx graph
	G = nx.Graph()
	# We fill in the vertices
	for i in range(n_cities):
		G.add_node(i)
	# We fill in all the edges
	for i in range(n_cities):
		for j in range(n_cities):
			if i<=j: 
				continue
			G.add_edge(i,j,weight=bias/weights[i][j])

	# Now we draw the graph
 	if pos==None:
 		pos = nx.layout.spring_layout(G)
	nx.draw_networkx_nodes(G,pos=pos,node_size=10)
	nx.draw_networkx_edges(G,pos,edgelist=edgelist,
                    width=edge_width,alpha=0.2,edge_color=edge_color)
	fig.canvas.draw()
	plt.show()
	return pos
tsp.plot = _plot_tsp

