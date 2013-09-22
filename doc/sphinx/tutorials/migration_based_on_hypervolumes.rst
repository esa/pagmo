.. _migration_based_on_hypervolumes:

================================================================
Migration based on hypervolume contribution
================================================================

In this tutorial we will cover some migration strategies that are based on the hypervolume computation.
There are in total 4 migration policies which are based on the hypervolume feature.

#. `PyGMO.migration.hv_greedy_s_policy`
#. `PyGMO.migration.hv_greedy_r_policy`
#. `PyGMO.migration.hv_best_s_policy`
#. `PyGMO.migration.hv_fair_r_policy`

For more information on migration policies, please visit the corresponding documentation page on :ref:`migration`.

Code below establishes an archipelago with the SMS-EMOA algorithm, using the `PyGMO.migration.hv_best_s_policy` and `PyGMO.migration.hv_fair_r_policy`.
As a comparison, the random migration policies `PyGMO.migration.random_s_policy` and `PyGMO.migration.random_r_policy` are also employed for the task.

.. code-block:: python

  from PyGMO import *

  def run_evolution(islands, prob):
    """
    Creates archipelago and proceeds with the evolution
    """

    # Create the archipelago, and push the islands
  	arch = archipelago(topology=topology.fully_connected())
  	for isl in islands:
  		arch.push_back(isl)
  
  	# Evolve for 130 steps
  	n_steps = 130
  	for s in xrange(n_steps):
  		print "Evolving archipelago, step %d/%d" % (s, n_steps)
  		arch.evolve(1)
  
  	# Merge all populations across the islands together
  	pop = population(prob)
  	for isl in arch:
  		for ind in isl.population:
  			pop.push_back(ind.cur_x)
  
  	print "Final P-Distance: ", prob.p_distance(pop)
  	prob.plot(pop)
  
  def main():
  	# Set up problem as DTLZ-3 with 3 objectives and the algorithm as SMS-EMOA
  	prob = problem.dtlz3(fdim=3)
  	alg = algorithm.sms_emoa(gen = 100)
  
  	# Construct the hv_best/fair migration policies
  	s_pol = migration.hv_best_s_policy(0.25, migration.rate_type.fractional)
  	r_pol = migration.hv_fair_r_policy(0.25, migration.rate_type.fractional)
  
  	# Construct the random policies
  	r_s_pol = migration.random_s_policy(0.25, migration.rate_type.fractional)
  	r_r_pol = migration.random_r_policy(0.25, migration.rate_type.fractional)
  
  	# Set up the archipelago
  	n_islands = 16
  	n_individuals = 64
  
  	# Create and evolve the archipelago using the hypervolume-based migration policies
  	isls_hv = [island(alg, prob, n_individuals, s_policy=s_pol, r_policy=r_pol) for i in xrange(n_islands)]
  	run_evolution(isls_hv, prob)
  
  	# Create and evolve the archipelago using the random migration policies
  	isls_rnd = [island(alg, prob, n_individuals, s_policy=r_s_pol, r_policy=r_r_pol) for i in xrange(n_islands)]
  	run_evolution(isls_rnd, prob)

  if __name__ == "__main__":
    main()

.. note::
 You can save the code above, and execute it by issuing the following in the command line: **python tutorial.py** (assuming the first argument is the name of the file).

After 130 evolutionary steps, the first scenario produces a population which has converged to a solution not far from the true pareto front.
The plot below is a result of the evolution of an archipelago using the hypervolume-based migration policies:

.. image:: ../images/tutorials/hv_best_fair_migration_policy.png
  :width: 750px

In case of the random migration policies, the individuals are still far from the optimal front, which suggests that the hypervolume-based migration policies might have helped in the establishing of the good solution.
Plot below is a result of the evolution of an archipelago using the random migration policies:

.. image:: ../images/tutorials/random_migration_policy.png
  :width: 750px

How does the migration work ?
=============================

We owe you an explanation on what had just happened behind the curtains of that archipelago migration.
The core idea of having the evolution happen on the archipelago are the occasional migrations - some individuals from given island
may occasionally be copied and sent to a neighbouring island.
The island to which the individuals have travelled is able to pick and choose the newly arrived immigrants, and use the information stored in their chromosome to advance the evolution further.

Hypervolume computation plays a significant role in establishing the *best* subset of individuals (these are the candidates for emigration), as well as the *worst* subset (which may be replaced by available set of immigrants). 

Hypervolume-based selection policy
----------------------------------

Let us assume a population of 10 individuals, and we want to select the set of immigrants of size 3.

.. image:: ../images/tutorials/hv_migration_selection.png
  :width: 750px

Hypervolume-based replacement policy
------------------------------------

.. image:: ../images/tutorials/hv_migration_replacement.png
  :width: 750px

