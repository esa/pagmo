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
import algorithm, migration, problem, topology, kep_toolbox
from copy import copy

def run_test():
	from PyGMO import problem, algorithm, island
	from numpy import mean, std
	number_of_islands = 5
	number_of_individuals = 30
	number_of_generations = 250

	prob_list = [problem.himmelblau(), problem.schwefel(10), problem.rastrigin(20), problem.griewank(5), problem.rosenbrock(4), problem.dejong(3), problem.michalewicz(5)]

	#Known solutions to the optimization problems
	prob_optimum = [0, 0, 0, 0, 0, 0, -4.687]

	algo_list = [algorithm.pso(number_of_generations), algorithm.de(number_of_generations,0.8,0.8),algorithm.sa_corana(number_of_generations,1,0.1), algorithm.ihs(number_of_generations), algorithm.cs(number_of_generations, 0.1, 0.25), algorithm.bee_colony(number_of_generations), algorithm.firefly(number_of_generations) ]

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
			print('                Best:\t' + str(min(best)[0] - prob_optimum[j]))
			print('                Mean:\t' + str(mean(best) - prob_optimum[j]))
			print('                Std:\t' + str(std(best)))

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


