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
__version__ = '1.0.0'

# For convenience, bring all core classes into the root namespace when importing *.
from core import *
__all__ += filter(lambda name: not name.startswith('_'),dir(core))

problem_list = problem._get_problem_list()
algorithm_list = algorithm._get_algorithm_list()
island_list = core._get_island_list()

# Fill up the __extensions__ variable with all detected extensions
__extensions__ = {'nlopt': False, 'gsl': False,'snopt': False,'ipopt': False,'gtop': False,'scipy': False,'networkx': False,'vpython': False, 'pykep': False}
if "nlopt" in str(algorithm._get_algorithm_list()):
	__extensions__['nlopt']=True
if "gsl" in str(algorithm._get_algorithm_list()):
	__extensions__['gsl']=True
if "snopt" in str(algorithm._get_algorithm_list()):
	__extensions__['snopt']=True
if "ipopt" in str(algorithm._get_algorithm_list()):
	__extensions__['ipopt']=True
if "cassini" in str(problem._get_problem_list()):
	__extensions__['gtop']=True
try:
	import scipy
	__extensions__['scipy']=True
except:
	pass
try:
	import networkx
	__extensions__['networkx']=True
except:
	pass
try:
	import visual
	__extensions__['vpython']=True
except:
	pass

try:
	import PyKEP
	__extensions__['pykep']=True
except:
	pass

def run_test(n_trials=200, pop_size = 20, n_gen = 500):
	from PyGMO import problem, algorithm, island
	from numpy import mean, std
	number_of_trials = n_trials
	number_of_individuals = pop_size
	number_of_generations = n_gen

	prob_list = [problem.schwefel(dim = 10), problem.rastrigin(dim = 10), problem.rosenbrock(dim = 10), problem.ackley(dim = 10), problem.griewank(dim = 10), problem.levy5(10)]
	if __extensions__['gtop']:
		prob_list.append(problem.cassini_1())
		prob_list.append(problem.gtoc_1())
		prob_list.append(problem.cassini_2())
		prob_list.append(problem.messenger_full())
		
	algo_list = [algorithm.pso(gen = number_of_generations), algorithm.de(gen = number_of_generations,xtol=1e-30, ftol=1e-30), algorithm.de_self_adaptive(gen = number_of_generations, restart=True,xtol=1e-30, ftol=1e-30), algorithm.de_1220(gen = number_of_generations, restart=True,xtol=1e-30, ftol=1e-30), algorithm.sa_corana(iter = number_of_generations*number_of_individuals,Ts = 1,Tf = 0.01), algorithm.ihs(iter = number_of_generations*number_of_individuals), algorithm.sga(gen = number_of_generations), algorithm.cmaes(gen = number_of_generations,xtol=1e-30, ftol=1e-30), algorithm.bee_colony(gen = number_of_generations/2)]
	print('\nTrials: ' + str(n_trials) + ' - Population size: ' + str(pop_size) + ' - Generations: ' + str(n_gen))
	for prob in prob_list:
		print('\nTesting problem: ' + prob.get_name() + ', Dimension: ' + str(prob.dimension) )
		print('With Population Size: ' +  str(pop_size) )
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

def example_1(n_trials=25, variant_adptv=2, restart=True):
	from PyGMO import problem, algorithm, island, archipelago
	from PyGMO.topology import fully_connected
	from numpy import mean, median
	results = list()
	prob = problem.messenger_full()
	de_variants = [11,13,15,17]
	algos = [algorithm.de_self_adaptive(gen=50,variant=v, restart=restart, variant_adptv=variant_adptv) for v in de_variants]
	
	for trial in range(n_trials):
		archi = archipelago(topology=fully_connected())
		for algo in algos:
			archi.push_back(island(algo,prob,25)) 
		print "Trial N: " + str(trial)
		archi.evolve(30)
		results.append(min([isl.population.champion.f[0] for isl in archi]))
	return (mean(results), median(results), min(results), max(results))
	
def example_2(algo=algorithm.de(1), prob = problem.rosenbrock(10), topo = topology.barabasi_albert(3,3), n_evolve = 100, n_isl = 1024, pop_size = 20, color_code='rank'):
	from PyGMO import problem, algorithm, island, archipelago
	from matplotlib.pyplot import savefig, close
	archi = archipelago(algo,prob,n_isl,pop_size,topology=topo)
	print "Drawing Initial Condition .. "
	pos = archi.draw(scale_by_degree=True,n_size=3,e_alpha=0.03, n_color = color_code)
	savefig('archi000', dpi = 72)
	close()
	for i in range(1,n_evolve):
		archi.evolve(1); 
		archi.join();
		print "Drawing"+ str(i) +  "-th evolution .. "
		pos = archi.draw(layout = pos, scale_by_degree=True,n_size=3,e_alpha=0.03, n_color = color_code)
		savefig('archi%03d' % i, dpi = 72);  
		close()

#def test_aco():
#	from PyGMO import problem, algorithm, island
#	from numpy import mean, std
#	number_of_islands = 5
#	number_of_individuals = 30
#	number_of_generations = 50
#	w = [ [0,1,100,1], [1,0,1,100], [100,1,0,1], [1, 100, 1, 0]]
#	prob_list = [problem.tsp(w)]
#	algo_list = [algorithm.aco(number_of_generations)]
#	for j in range(0,len(prob_list)):
#		print('Testing problem: ' + str(type(prob_list[j])) + ', Dimension: ' + str(prob_list[j].dimension))
#		for algo in algo_list:
#			print('        Testing algorithm: ' + str(algo))
#			best = []
#			best_x = []
#			for i in range(0,number_of_islands):
#				isl = island(prob_list[j],algo,number_of_individuals)
#				best.append(isl.population.champion.f)
#				best_x.append(isl.population.champion.x)
#				print('                Best fitness:\t' + str(best[i]))
#				print('                Best solution:\t' + str(best_x[i]))

