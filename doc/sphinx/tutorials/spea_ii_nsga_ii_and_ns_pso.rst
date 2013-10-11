.. _spea_ii_nsga_ii_and_ns_pso:

================================================================
Using NSGA-II, SPEA2 and NS-PSO
================================================================

We will now introduce 3 more multi-objective optimization algorithms.

Let's start with **NSGA-II**. NSGA-II is a non-dominated sorting based multi-objective evolutionary algorithm. It generates offspring with crossover and mutation and select the next generation according to non-dominated sorting and crowding distance comparison.
As for PADE it is possible to solve a multi-objective optimization problem with NSGA-II as follow

.. code-block:: python

	In [1]: from PyGMO import *
	In [2]: prob = problem.zdt(1)
	In [3]: alg = algorithm.nsga_II()
	In [4]: pop = population(prob, 100)
	In [5]: pop = alg.evolve(pop)
	In [6]: pop.plot_pareto_fronts()

It is possible to specify the number of generations to run the algorithm for, crossover and mutation probability as well as the mutation and crossover distribution index, as follow

.. code-block:: python

	In [7] alg = algorithm.nsga_II(gen = 100, cr = 0.95, eta_c = 10, m = 0.01, eta_m = 50)

This will run the NSGA-II algorithms for 100 generations, with a crossover probability of 0.95, mutation probability of 0.01, crossover distribution index of 10 and mutation distribution index of 50.

We now introduce **NSPSO**. Non-dominated Sorting Particle Swarm Optimizer (NSPSO) is a modified version of PSO for multi-objective optimization. It extends the basic ideas of PSO by making a better use of personal bests and offspring for non-dominated comparison. As for PSO it is possible to set:

 * *C1* and *C2* (the magnitude of the force to apply towards respectively the personal best and the global best of a particle)
 * *CHI* the velocity scaling factor
 * *m_v_coeff* the velocity coefficient (determining the maximum allowed particle velocity)
 * *minW* and *maxW* which defines in which range the inertia weight will be adapted throughout the run.

Here we see how those parameters can be set when instantiating the algorithm

.. code-block:: python

	In [8]: alg = algorithm.nspso(gen = 10, minW = 0.4, maxW = 1.0, C1 = 2.0, C2 = 2.0, CHI = 1.0, v_coeff = 0.5)

NSPSO selects the global best for each particles among non-dominated particles. The non-dominated particles are sorted according to one niching method (crowding distance, niche count or maxmin) and the leader is selected among the best ones. The parameter *leader_selection_range* define which fraction of the non-dominated particles to consider for selection as global best. 

In the following code leader_selection_range is set to 20. That means that the global best for each particle is randomly selected among the top 20% (according to the niching method) non-dominated individuals.

.. code-block:: python
	
	In [9]: alg = algorithm.nspso(leader_selection_range = 20)

It is possible to choose between 3 different niching methods: CROWDING DISTANCE, NICHE COUNT and MAXMIN. This can be done as follow:

.. code-block:: python

	In [10]: alg = algorithm.nspso(diversity_mechanism = algorithm.nspso.CROWDING_DISTANCE)
	In [11]: alg = algorithm.nspso(diversity_mechanism = algorithm.nspso.NICHE_COUNT)
	In [12]: alg = algorithm.nspso(diversity_mechanism = algorithm.nspso.MAXMIN)

The last multi-objective optimization algorithm we introduce is **SPEA2**. In the Strength Pareto Evolutionary Algorithm (SPEA2) the quality of an individual is measured taking into consideration its pareto strength and its distance to its K-th neighbour, where *K = sqrt(pop size + archive size)*.
It uses the same mutation and crossover operators of NSGA-II, so as for the former it is possible to specify the number of generations to run the algorithm for, crossover and mutation probability as well as the mutation and crossover distribution index, as follow.

.. code-block:: python

	In [13]: alg = algorithm.spea2(gen = 100, cr = 0.95, eta_c = 10, m = 0.01, eta_m = 50)

SPEA2 uses an external archive in which are stored the non dominated solutions found so far. The size of the archive is kept constant throughout the run by mean of a truncation operator taking into consideration the distance of each individual to its closest neighbours. 

It is possible to set the archive size in two ways. If *archive_size* is set to 0 (which is the default behaviour), then archive size is adapted to be equal to the population which is evolving. In the following case, for example, the archive size will be set equal to 100.

.. code-block:: python

	In [13]: alg = algorithm.spea2(archive_size = 0)
	In [14]: pop = population(prob, 100)
	In [15]: pop = alg.evolve(pop)

Otherwise it is possible to set the archive size to any size, up to the population size.

.. code-block:: python

	In [16]: alg = algorithm.spea2(archive_size = 20)
