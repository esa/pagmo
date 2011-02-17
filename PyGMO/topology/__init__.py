# -*- coding: utf-8 -*-
from _topology import *

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
			retval.add_edges_from([(i,n) for n in self.get_adjacent_vertices(i)])
		else:
			retval.add_node(i)
	return retval

def _draw(self):
	"""
	Draw topology using the draw() command from networkx.
	"""
	try:
		import networkx as nx
	except ImportError:
		raise ImportError('Could not import the networkx module.')
	if self.number_of_vertices == 0 or self.number_of_vertices == 1:
		raise ValueError('Cannot draw topology with one single vertex or less.')
	nx.draw(self.to_networkx())

_topology._base.to_networkx = _to_networkx
_topology._base.draw = _draw
