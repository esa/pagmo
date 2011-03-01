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

import core, algorithm, migration, problem, topology, test

__doc__ = 'PyGMO is a pretty cool guy. eh kills aleins and doesnt afraid of anything...'
__all__ = ['core', 'algorithm', 'migration', 'problem', 'topology', 'test']

# For convenience, bring all core classes into the root namespace when importing *.
from core import *
__all__ += filter(lambda name: not name.startswith('_'),dir(core))

problem_list = problem._get_problem_list()
algorithm_list = algorithm._get_algorithm_list()
island_list = core._get_island_list()

def run_test():
	from PyGMO import problem, algorithm, island
	from numpy import mean, std
	number_of_trials = 200
	number_of_individuals = 20
	number_of_generations = 500

	prob_list = [problem.schwefel(10), problem.rastrigin(10), problem.rosenbrock(10), problem.ackley(10), problem.griewank(10)]
	algo_list = [algorithm.pso(number_of_generations), algorithm.de(number_of_generations,0.8,0.8,2),algorithm.sa_corana(number_of_generations*number_of_individuals,1,0.1), algorithm.ihs(number_of_generations*number_of_individuals), algorithm.sga(number_of_generations,0.8,0.1)]

	for prob in prob_list:
		print('\nTesting problem: ' + str(type(prob)) + ', Dimension: ' + str(prob.dimension))
		for algo in algo_list:
			print(' ' + str(algo))
			best = []
			best_x = []
			for i in range(0,number_of_trials):
				isl = island(algo,prob,number_of_individuals)
				isl.evolve(1)
				isl.join()
				best.append(isl.population.champion.f)
				best_x.append(isl.population.champion.x)
			print(' Best:\t' + str(min(best)[0]))
			print(' Mean:\t' + str(mean(best)))
			print(' Std:\t' + str(std(best)))

def test_aco():
	from PyGMO import problem, algorithm, island
	from numpy import mean, std
	number_of_islands = 5
	number_of_individuals = 30
	number_of_generations = 50
	w = [ [0,1,100,1], [1,0,1,100], [100,1,0,1], [1, 100, 1, 0]]
	prob_list = [problem.tsp(w)]
	algo_list = [algorithm.aco(number_of_generations)]
	for j in range(0,len(prob_list)):
		print('Testing problem: ' + str(type(prob_list[j])) + ', Dimension: ' + str(prob_list[j].dimension))
		for algo in algo_list:
			print('        Testing algorithm: ' + str(algo))
			best = []
			best_x = []
			for i in range(0,number_of_islands):
				isl = island(prob_list[j],algo,number_of_individuals)
				isl.evolve(1)
				isl.join()
				best.append(isl.population.champion.f)
				best_x.append(isl.population.champion.x)
				print('                Best fitness:\t' + str(best[i]))
				print('                Best solution:\t' + str(best_x[i]))


