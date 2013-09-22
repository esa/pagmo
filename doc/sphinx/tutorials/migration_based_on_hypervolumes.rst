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

Let us assume an island with 10 individuals, out of which we want to determine a set of 4 emigrants (outgoing individuals).
First step would be computing the contributions of each individual, and then selecting a subset of 4 individuals that contributed the most volume.
Plot on the left visualizes the computed exclusive contributions of 10 individuals.
Plot on the right, are the same individuals ordered descending by their contribution.
First four individuals in new ordering are candidates for emigrants.

.. image:: ../images/tutorials/hv_migration_selection.png
  :width: 750px

Description above is almost exactly what happens in `PyGMO.migration.hv_best_s_policy`.
The main difference in what hypervolume-based "best" policy does is not computing the contributions for the whole population in a single call, but doing it *per-front*, starting from the first one.
This is mainly a precaution for selecting the best emigrants possible in the early stages of the algorithm.
Since many points are dominated (thus having a zero contribution) we would act more or less arbitrarily if we were to just rely on that measure.

In order to do that we first compute the contributions among the individuals in the first front (as every individual in other front has a contribution of 0).
If the requested number of individuals for emigration was not met yet, we remove the first front from the temporary population, and continue to fill up the list of emigrants with the greatest contributors from the original second front. This continues until we find the requested number of emigrants.

Although the general idea of `PyGMO.migration.hv_greedy_s_policy` is the same, there are few differences among them:
We don't establish the ordering of individuals also according to the exclusive hypervolume, but this time this is done iteratively.
Instead of computing the contributions of all individuals at once (see `PyGMO.hypervolume.contributions`), we explicitly request for the computation of the greatest contributor (see `PyGMO.hypervolume.greatest_contributor`).
The *per-front* policy also applies here.
After the greatest contributor is found we immediately remove it from the temporary population, and continue computation of the next greatest contributor (this continues until we have found sufficient number of greatest contributors).
It's easy to notice that once given point is removed from the population, the exclusive contributions of other points may change.

Hypervolume-based replacement policy
------------------------------------

Replacement policy works in a similar fashion, except this time it's the *K* least contributors of the merged sets of immigrants (incoming individuals) and islanders (individuals belonging to given island).
Plot on the left visualizes a set of 10 islanders and 5 immigrants merged together.
In thus created set of 15 individuals, we determine a subset of 5 least contributors, either by computing all contributions as once (`PyGMO.migration.hv_fair_r_policy`) or iteratively (`PyGMO.migration.hv_greedy_r_policy`) by removing each least contributor (`PyGMO.hypervolume.least_contributor`) once it was established.
Plot on the right visualizes the ordered set, out of which 5 least contributors were selected.
Since there are 3 islanders in the set of 5 least contributors, it is possible to make 3 fair replacements: 3 *discarded* islanders (crossed-over bar) and 3 *non-discarded* immigrants.

.. image:: ../images/tutorials/hv_migration_replacement.png
  :width: 750px

.. note::
 The *per-front* policy also applies here. Least contributors are established first from the **last** front of the population, progressing upwards to the individuals in the first front.

Additional thing both greedy and fair migrations do is filtering-out the duplicated individuals.
It is likely that strong individuals will start to prevail among the immigrants.
If given individual is already on given island, we would like to make sure it is not added to the population, as having a diverse chromosome pool is crucial.
In order to make sure, no scenario like that happens, the set of immigrants is matched against a set of islanders. Any duplicated immigrants are thus discarded up-front.
