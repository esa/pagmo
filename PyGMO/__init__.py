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

from _PyGMO import *
from copy import copy

class __vector_init:
	def build(self, iterable, t):
		if iterable == None:
			return
		if not getattr(iterable, '__iter__', False):
			raise TypeError, 'I need an iterable object for initialisation.'
		for i in iterable: self.append(t(i))


class vector_double(_PyGMO.__base_vector_double,__vector_init):
	def __init__(self, iterable = None):
		super(type(self), self).__init__()
		self.build(iterable,float)

class vector_size_t(_PyGMO.__base_vector_size_t,__vector_init):
	def __init__(self, iterable = None):
		super(type(self), self).__init__()
		self.build(iterable,int)

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
