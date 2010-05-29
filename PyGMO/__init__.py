# -*- coding: iso-8859-1 -*-
# Copyright (C) 2004-2009 The PaGMO development team,
# Advanced Concepts Team (ACT), European Space Agency (ESA)
# http://apps.sourceforge.net/mediawiki/pagmo
# http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers
# http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits
# act@esa.int
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
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
import algorithm, migration, problem, topology
from copy import copy

def run_test():
	from PyGMO import problem, algorithm, island
	from numpy import mean, std
	prob_list = [problem.schwefel(10),problem.rastrigin(10),problem.rosenbrock(10), problem.ackley(10)]
	algo_list = [algorithm.pso(500), algorithm.de(500,0.8,0.8),algorithm.sa_corana(10000,1,0.1), algorithm.ihs(500), algorithm.gsl_nm(500,0.5,0.5), algorithm.sga(500, 0.8, 0.1), algorithm.cs(500, 0.1, 0.25) ]
	for prob in prob_list:
		print('Testing problem: ' + str(type(prob)) + ', Dimension: ' + str(prob.dimension))
		for algo in algo_list:
			print('        Testing algorithm: ' + str(algo))
			best = []
			for i in range(0,100):
				isl = island(prob,algo,20)
				isl.evolve(1)
				isl.join()
				best.append(isl.population.champion.f)
			print('                Best:\t' + str(min(best)))
			print('                Mean:\t' + str(mean(best)))
			print('                Std:\t' + str(std(best)))
