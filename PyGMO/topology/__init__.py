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

def _draw(self, layout = "spring", n_size = 15, scale_by_degree = False, n_color = 'blue', n_alpha = 0.5, e_alpha = 0.1):
	"""
	Draw topology using the draw() command from networkx.

        USAGE: topology.draw(layout = "spring", n_size = 15, scale_by_degree = False, n_color = 'blue', n_alpha = 0.5, e_alpha = 0.1)

        * layout: Network layout. Can be 'spring' or 'circular'.
        * n_size: The size of nodes. Becomes scaling factor when scale_by_degree=True.
        * scale_by_degree: When True, nodes will be sized proportional to their degree.
        * n_color: Node color.
        * n_alpha: Transparency of nodes. Takes value between 0 and 1.
        * e_elpha: Transparency of edges. Takes value between 0 and 1.
	"""
	try:
		import networkx as nx
	except ImportError:
		raise ImportError('Could not import the networkx module.')
        try:
                import matplotlib as mpl
        except ImportError:
                raise ImportError('Could not improt the MatPlotLib module.')
	if self.number_of_vertices == 0 or self.number_of_vertices == 1:
		raise ValueError('Cannot draw topology with one single vertex or less.')

        G = self.to_networkx()
        node_sizes = range(nx.number_of_nodes(G))
        for i in range(nx.number_of_nodes(G)):
                if scale_by_degree:
                        node_sizes[i] = nx.degree(G,i)*n_size
                else:
                        node_sizes[i] = n_size

        if layout == "spring":
                pos = nx.spring_layout(self.to_networkx())
        if layout == "circular":
                pos = nx.circular_layout(self.to_networkx())

        mpl.pyplot.figure()
        nx.draw_networkx_edges(self.to_networkx(),pos,alpha=e_alpha,arrows=False)
        nx.draw_networkx_nodes(self.to_networkx(),pos,node_size=node_sizes,node_color=n_color,alpha=n_alpha,with_labels=False)
        mpl.pyplot.axis('off')
        mpl.pyplot.show()


def _draw_degree_distribution(self, style = '.'):
        """
        Plot the degree dstribution on a loglog scale

        USAGE topology.draw_degree_distribution(style = 'r.')

        style: MatPlotLib line style.
        """
        try:
                import matplotlib as mpl
        except ImportError:
                raise ImportError('Could not improt the MatPlotLib module.')
        dd = self.get_degree_distribution();
        mpl.pyplot.figure()
        mpl.pyplot.loglog(dd,style)
        mpl.pyplot.xlabel('k - number of links')
        mpl.pyplot.show('p(k) - probability of k')

_topology._base.to_networkx = _to_networkx
_topology._base.draw = _draw
_topology._base.draw_degree_distribution = _draw_degree_distribution
