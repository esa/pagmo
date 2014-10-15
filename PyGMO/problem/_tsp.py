from PyGMO.problem._problem import tsp,_tsp_encoding


# Renaming and placing the enums
tsp.encoding = _tsp_encoding

def _tsp_ctor(self, matrix=[[0, 1, 2], [1, 0, 5], [2, 5, 0]], type="full"):
    """
    Constructs a Travelling Salesman problem (Constrained Integer Single-Objective)

    USAGE: problem.tsp(matrix = [0,1,2],[1,0,5],[2,5,0])

    * matrix: inter-city distances (symmetric matrix)
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(matrix)
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
	edgelist = self.chromosome2cities(x)
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

