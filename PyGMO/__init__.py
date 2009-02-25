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

def adaptive_optimization(prob,topo,algo_list,pop_size,retval):
	from PyGMO import archipelago,island
	from copy import copy
	from numpy import array
	for iter in range(0,200):
		# Build the initial archipelago: same island, different algorithms.
		a = archipelago(prob,topo)
		isl = island(prob,algo_list[0],pop_size)
		for i in algo_list:
			a.append(isl)
			a[-1].algorithm = i
		counter = 0;
		while True:
			if counter > 20:
				break
			# Calculate the current best individuals for each island.
			cur_best = array([i.best().fitness for i in a])
			# Evolve with different algorithm.
			a.evolve_t(12000)
			a.join()
			# New best individuals.
			new_best = array([i.best().fitness for i in a])
			# Calculate the variation of best individuals after the evolution.
			delta = cur_best - new_best
			# The best algorithm is the one that brought the greater improvement.
			if max(delta) == 0.:
				print "No difference in delta, continuing with different algos."
				counter += 1
				continue
			else:
				counter = 0
			idx = list(delta).index(max(delta))
			print "Best algorithm = ", algo_list[idx].__repr__()
			print delta
			# Copy the best island all over the archipelago.
			for i in range(0,len(a)): a[i] = a[idx]
			# Evolve the archipelago.
			a.evolve_t(8000)
			a.join()
			# Find the best individuals.
			best = [i.best().fitness for i in a]
			# Copy the best island all over the archipelago and change the algorithms.
			idx = best.index(min(best))
			for i in range(0,len(a)):
				a[i] = a[idx]
				a[i].algorithm = algo_list[i]
			print a[0].best().fitness
		retval.append(a)
		print "New iteration."

def vector(x, t = None):
	import PyGMO.core as core
	import re
	l = dir(core)
	p = re.compile('vector_.*')
	vector_types = []
	if t == None:
		for i in l:
			if re.match(p,i):
				vector_types.append(getattr(core, i))
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

def __arch_prune(arch,perc,display = False):
	if display:
		from matplotlib.pylab import subplot,plot,xlim
	from math import sqrt
	from PyGMO import vector_double
	if not isinstance(perc,int) or perc < 0 or perc > 100:
		raise ValueError('percentile must be an integer in the [0,100] range')
	ind_list = [isl.best() for isl in arch]
	ind_list = sorted(ind_list, key = lambda i: i.fitness)[0:(len(ind_list)*perc)/100]
	if len(ind_list) == 0:
		raise ValueError('the given percentile results in an empty list of best individuals')
	prob = arch.problem
	p_dimension = prob.dimension
	if display:
		edge_size = int(sqrt(p_dimension))
		height = edge_size
		width = edge_size
		if height * width != p_dimension:
			height += 1
	new_bounds = (vector_double(),vector_double())
	for i in range(0,p_dimension):
		if display:
			subplot(width,height,i+1)
			plot([ind.decision_vector[i] for ind in ind_list], [ind.fitness for ind in ind_list],'o')
			xlim((prob.lb[i],prob.ub[i]))
		old_width = prob.ub[i] - prob.lb[i]
		new_lb = min([ind.decision_vector[i] for ind in ind_list]) - old_width * .1
		new_ub = max([ind.decision_vector[i] for ind in ind_list]) + old_width * .1
		new_bounds[0].append(max(prob.lb[i],new_lb))
		new_bounds[1].append(min(prob.ub[i],new_ub))
	retval = arch.problem
	retval.lb = new_bounds[0]
	retval.ub = new_bounds[1]
	return retval

archipelago.make_neato = __arch_make_neato
archipelago.prune = __arch_prune
