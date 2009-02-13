# Copyright (C) 2007, 2008 by Francesco Biscani
# bluescarni@gmail.com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

from core import *
import algorithm, problem, topology
from copy import copy

def vector(x, t = None):
	import _PyGMO as PyGMO
	import re
	l = dir(PyGMO)
	p = re.compile('vector_.*')
	vector_types = []
	if t == None:
		for i in l:
			if re.match(p,i):
				vector_types.append(getattr(PyGMO, i))
		if len(vector_types) == 0:
			raise TypeError('No vector classes in PyGMO.')
	else:
		vector_types.append(t)
	retval = None
	for i in vector_types:
		retval = i()
		if getattr(x, '__iter__', False):
			try:
				for j in x:
					retval.append(j)
				return retval
			except TypeError:
				pass
		else:
			try:
				retval.append(x)
				return retval
			except TypeError:
				pass
	raise TypeError("No suitable vector class found in PyGMO.")

def __arch_make_neato(arch,directed = True):
	if directed:
		graph_t = 'digraph'
		graph_sep = '->'
	else:
		graph_t = 'graph'
		graph_sep = '--'
	t_str = arch.topology.__repr__()
	t_dict = {}
	for i in t_str.split():
		tmp = i.split('->')
		t_dict[int(tmp[0])] = [int(j) for j in tmp[1].split(',')]
	max_n_edges = max([len(t_dict[i]) for i in t_dict])
	min_n_edges = min([len(t_dict[i]) for i in t_dict])
	print 'Max number of edges = ', max_n_edges
	print 'Min number of edges = ', min_n_edges
	retval = 'strict ' + graph_t + ' foo {\n'
	retval += 'edge [len=' + str(max_n_edges) + '];\n'
	for i in t_dict:
		retval += str(i) + ' [height=' + str(len(t_dict[i])/5.) + ',width=' + str(len(t_dict[i])/5.) + ',label=' + str(len(t_dict[i])) + ',fontsize=' + \
			str(len(t_dict[i])*10) + '];\n'
		#retval += 'edge [w=' + str(len(t_dict[i])*10./max_n_edges) + '];\n'
		for j in t_dict[i]:
			retval += str(i) + graph_sep + str(j) + ';\n'
	retval += '}'
	return retval

archipelago.make_neato = __arch_make_neato
