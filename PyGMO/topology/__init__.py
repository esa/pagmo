# -*- coding: iso-8859-1 -*-
from _topology import *

def _to_networkx(self):
	try:
		import networkx as nx
	except ImportError:
		raise ImportError('Could not import the networkx module.')
	retval = nx.DiGraph()
	for i in range(self.number_of_vertices):
		if self.get_num_adjacent_vertices(i):
			retval.add_edges_from([(i,n) for n in self.get_adjacent_vertices(i)])
		else:
			retval.add_edge((i,None))
	return retval

def _draw(self):
	try:
		import networkx as nx
	except ImportError:
		raise ImportError('Could not import the networkx module.')
	nx.draw(self.to_networkx())

_topology._base.to_networkx = _to_networkx
_topology._base.draw = _draw
