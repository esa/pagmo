Testing Algorithms
==================

PaGMO can be used to test the performance of algorithms on problem sets. 
We show how simple is to do so using a basic set-up..

.. code-block:: python

   from PyGMO import *
   trials = 50
   individuals = 20
   generations = 500
   prob_list = [problem.schwefel(dim = 10), problem.rastrigin(dim = 10), problem.rosenbrock(dim = 10)]
   algo_list = [algorithm.pso(gen = generations), algorithm.de(gen = generations)]

We have instantiated the problems and algorithms we want to test. We have used default parameters
and variants implicitly, but these could be specified as well and compared. We now run the tests ... 
note that we use the island rather than the archipelago, thus the test WILL NOT run in multiple threads.
One could speed things up (in a multiple CPU environment) byt using the archipelago
and an unconnected topology.

.. code-block:: python

   for prob in prob_list:
	for algo in algo_list:
		best = []
		best_x = []
		for i in range(0,trials):
			# Here we create the optimization problem .....
			isl = island(algo,prob,individuals)
			# And here we solve it
			isl.evolve(1)
			isl.join()
			best.append(isl.population.champion.f)
			best_x.append(isl.population.champion.x)
	#Here you can store the results ....

The final script (included in PyGMO) is:

.. code-block:: python

   def run_test(n_trials=200, pop_size = 20, n_gen = 500):
	from PyGMO import problem, algorithm, island
	from numpy import mean, std
	number_of_trials = n_trials
	number_of_individuals = pop_size
	number_of_generations = n_gen

	prob_list = [problem.schwefel(dim = 10), problem.michalewicz(dim = 10), problem.rastrigin(dim = 10), problem.rosenbrock(dim = 10), problem.ackley(dim = 10), problem.griewank(dim = 10)]
	if __extensions__['gtop']:
		prob_list.append(problem.cassini_1())
		prob_list.append(problem.cassini_2())
		prob_list.append(problem.gtoc_1())
		prob_list.append(problem.rosetta())
		prob_list.append(problem.messenger_full())
		prob_list.append(problem.tandem(prob_id = 6, max_tof = 10))
		
	algo_list = [algorithm.pso(gen = number_of_generations), algorithm.de(gen = number_of_generations,xtol=1e-30, ftol=1e-30), algorithm.de_self_adaptive(gen = number_of_generations, restart=True, variant_adptv=2,xtol=1e-30, ftol=1e-30), algorithm.de_1220(gen = number_of_generations, restart=True, variant_adptv=2,xtol=1e-30, ftol=1e-30), algorithm.sa_corana(iter = number_of_generations*number_of_individuals,Ts = 1,Tf = 0.01), algorithm.ihs(iter = number_of_generations*number_of_individuals), algorithm.sga(gen = number_of_generations), algorithm.cmaes(gen = number_of_generations,xtol=1e-30, ftol=1e-30), algorithm.bee_colony(gen = number_of_generations/2)]
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

You may check the output created by PyGMO.run_test() online at `PaGMO Tests <http://sourceforge.net/apps/mediawiki/pagmo/index.php?title=Tests>`_
