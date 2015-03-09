.. _migration:

Migration 
=========

The algorithm
-------------

Migration, in PyGMO, happens asynchronously in each :class:`PyGMO.island` between calls of the evolve() method 
of the  :class:`PyGMO.archipelago` where the :class:`PyGMO.island` has been pushed back. Each islands maintains a 
database of outgoing (or incoming) migrants. The various databases are used before evolution to replace
some :class:`PyGMO.individual` in the :class:`PyGMO.island` and are then updated after evolution with the
new selected migrants.

The algorithm is rather complex and the user does not need to know/understand its details as PyGMO sets defaults
values for all of its many parameters. These are:

Migration Rates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Default value: 1

These parameters defines the number of individuals that are selected from each :class:`PyGMO.island` for *migration* as well as
the number of migrants that will be considered for insertion in each :class:`PyGMO.island`.
To Migration Rates are set by the 's_policy' and 'r_policy' kwarg in the :class:`PyGMO.island` constructor. This can be done by specifying
the absolute number of individuals (migration.rate_type.absolute) or the fraction of the :class:`PyGMO.population` individuals
(migration.rate_type.fractional)

.. code-block:: python

   from PyGMO import *
   prob = problem.schwefel(15)
   algo = algorithm.de(gen = 100) #instantiates differential evolution with default params and 100 generations
   selection = migration.best_s_policy(0.10,migration.rate_type.fractional)
   replacement = migration.fair_r_policy(0.25,migration.rate_type.fractional)
   isl = island(algo,prob,s_policy = selection, r_policy = replacement)

Migration Direction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Default value: migration_direction.destination

In PyGMO the asynchronous *migration* is implemented by keeping a migrants database on each island. Then one of the following
options can be followed:

* 'destination': The internal migrants database stores, for each :class:`PyGMO.island`, the individuals that are meant to migrate from that island. Before each evolution, the :class:`PyGMO.island` will get migrating individuals from those made available by the islands connecting to it. After each evolution, the :class:`PyGMO.island` will update its list of best individuals in the database.

* 'source': The internal migrants database stores for each :class:`PyGMO.island` the individuals that are
   meant to migrate to that :class:`PyGMO.island`. Before each evolution, an :class:`PyGMO.island` will
   check if individuals destined to it are available in the database, and, in such case will, migrate over
   incoming individuals before starting evolution. After each evolution, the :class:`PyGMO.island` will
   place its candidate individuals for emigration in the database slots of the island(s) to which it connects.

The *migration* direction is set by the 'migration_direction' kwarg in the :class:`PyGMO.archipelago` constructor

.. code-block:: python

   from PyGMO import *
   prob = problem.schwefel(15)   
   algo = algorithm.de(gen = 100) #instantiates differential evolution with default params and 100 generations
   direction = migration_direction.source
   archi = archipelago(migration_direction = direction)


.. _distribution_type_label:

Migration Distribution Type 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Default value: distribution_type.point_to_point

When *migration* happens one has to decide which of the connected islands contributes to the event. This is decided
by the distribution type that can be one of the following:

* 'point to point': only one of the neighbourghing islands, selected at random, is sending (or receiving) the individuals

* 'broadcast': all neighbourghing islands are sending (or receiving) the individuals

The migration distribution type is set by the 'distribution_type' kwarg in the :class:`PyGMO.archipelago` constructor

.. code-block:: python

   from PyGMO import *
   prob = problem.schwefel(15)
   algo = algorithm.de(gen = 100) #instantiates differential evolution with default params and 100 generations
   distribution = distribution_type.broadcast
   archi = archipelago(distribution_type = distribution)

Migration Selection Policy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Default value: migration.best_s_policy(1)

The selection policy is the object responsible to choose out of a :class:`PyGMO.population` the individuals that will migrate. All selection policies derive from the same base class and currently a few are implemented:

* 'migration.best_s_policy': Selects the best individuals for a single-objective optimization problem. For a multi-objective optimization problem, an individual is considered better as another individual if it has a lower non-domination rank or - if the non-domination ranks of both individuals are equal - a higher crowding distance.

* 'migration.best_kill_s_policy': The same as 'migration.best_s_policy' but every selected individual gets reinitialized in the originating population.

* 'migration.random_s_policy': Individuals are selected uniformly at random.

* 'migration.hv_greedy_s_policy': Select the best individuals for a single-objective optimization problem. For a multi-objective optimization problem, and individual is considered better than another individual if its exclusive contribution to the hypervolume (see `PyGMO.hypervolume` for more details) is greater. The set of best individuals is created iteratively - after each selection of the individual, it is removed from the population so it does not diminish the contributions of other individuals.

* 'migration.hv_best_s_policy': Select the best individuals for a single-objective optimization problem. For a multi-objective optimization problem, and individual is considered better than another individual if its exclusive contribution to the hypervolume (see `PyGMO.hypervolume` for more details) is greater. The main distinction between this policy and `PyGMO.hv_best_s_policy` is computing all the contributions at once (without the removal step).


The selection policy is set by the 's_policy' kwarg in the :class:`PyGMO.island` constructor

.. code-block:: python

   from PyGMO import *
   prob = problem.schwefel(15)
   algo = algorithm.de(gen = 100) #instantiates differential evolution with default params and 100 generations
   best2 = migration.best_s_policy(2) # two individuals will be selected as the best
   isl = island(algo,prob,s_policy = best2)

Migration Replacement Policy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Default value: migration.fair_r_policy(1)

The replacement policy is the object responsible to substitute the individuals in a population with the
migrants. All replacement policies derive from the same base class and currently a few are implemented:

* 'migration.fair_r_policy': simply replaces the worst individuals in the island  with the best of the incoming migrants. This is subject to the added condition that the migrants are better. For multi-objective optimization problems, an individual is considered better than another individual if it has a lower non-domination rank or - if the non-domination ranks are equal - a higher crowding distance (compare with 'migration.best_s_policy')

* 'migration.random_r_policy': replaces random individuals in the island with random incoming migrants

* 'migration.worst_r_policy': replaces the worst individuals in the island with the best of the incoming migrants. In a multi-objective setting, the meaning of *better* is like in 'migration.fair_r_policy' or 'migration.best_s_policy').

* 'migration.hv_greedy_r_policy': Replaces a the worst individuals in the island with the best of the incoming immigrants.
  The distinction between individuals is made based on their exclusive contribution to the hypervolume (see `PyGMO.hypervolume` for more details).
  Both sets are determined iteratively - set of worst islanders is determined by choosing the least contributor among them, and then removing it from the population in order to prevent it from diminishing the contributions of other individuals.
  Likewise, the set of best immigrants is determined by their exclusive contribution to the hypervolume in an iterative fashion, except this time the greatest contributor is chosen.

* 'migration.hv_fair_r_policy': 
  The distinction between individuals is made based on their exclusive contribution to the hypervolume (see `PyGMO.hypervolume` for more details).
  Both sets are determined by computing the contributions to the hypervolume at once, without the removal step (as opposed to the 'migration.hv_greedy_r_policy').

The replacement policy is set by the 'r_policy' kwarg in the island constructor

.. code-block:: python

   from PyGMO import *
   prob = problem.schwefel(15)
   algo = algorithm.de(gen = 100) #instantiates differential evolution with default params and 100 generations
   random4 = migration.random_policy(4) # four individuals will be selected at random 
		 		        # from the migrants and will replace random individuals
   isl = island(algo,prob,s_policy = best2)


Migration Topology
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Default value: migration.unconnected()

The migration topology determines which island will be connected to which island. It also takes care that when an island is pushed back into an
archipelago, the topological properties of the resulting new connectivity graph are left unchanged.
It is set by the 'topology' kwarg in the archipelago constructor

.. code-block:: python

   from PyGMO import *
   prob = problem.schwefel(15)
   algo = algorithm.de(gen = 100) #instantiates differential evolution with default params and 100 generations
   topo = topology.ring()
   archi = archipelago(algo,prob,topology = topo)

The Classes
---------------------

.. class:: PyGMO.migration.best_s_policy([n=1, type = migration.rate_type.absolute])

   A selection policy that selects the n best :class:`PyGMO.individual` in
   the :class:`PyGMO.island`'s :class:`PyGMO.population`. If type is migration.rate_type.fractional then n, in [0,1], is interpreted
   as the fraction of the population to be selected. This class is used exclusively in the :class:`PyGMO.island` 
   constructor as a possible kwarg for the key 's_policy'

   .. code-block:: python

      from PyGMO import *
      prob = problem.griewank(5)
      algo = algorithm.bee_colony(gen = 10) #instantiates artificial bee colony with default params and 10 generations
      best2 = migration.best_s_policy(2)
      best50pc = migration.best_s_policy(0.5,migration.rate_type.fractional)
      isl1 = island(algo,prob,10,s_policy = best2)  #2 of the best individuals will migrate
      isl2 = island(algo,prob,32,s_policy = best50pc) #50% of 32 (i.e. 16) best individuals will migrate

.. class:: PyGMO.migration.best_kill_s_policy([n=1, type = migration.rate_type.absolute])

   A selection policy that selects the n best :class:`PyGMO.individual` in
   the :class:`PyGMO.island`'s :class:`PyGMO.population` and kills them in the original population so that
   only the migrant will survive. A new random individual will replace the migrant in the original population
   If type is migration.rate_type.fractional then n, in [0,1], is interpreted
   as the fraction of the population to be selected. This class is used exclusively in the :class:`PyGMO.island` 
   constructor as a possible kwarg for the key 's_policy'

   .. code-block:: python

      from PyGMO import *
      prob = problem.griewank(5)
      algo = algorithm.bee_colony(gen = 10) #instantiates artificial bee colony with default params and 10 generations
      best2 = migration.best_kill_s_policy(2)
      best50pc = migration.best_s_policy(0.5,migration.rate_type.fractional)
      isl1 = island(algo,prob,10,s_policy = best2)  #2 of the best individuals will migrate and be reinitialized in pop
      isl2 = island(algo,prob,32,s_policy = best50pc) #50% of 32 (i.e. 16) best individuals will migrate

.. class:: PyGMO.migration.random_s_policy([n=1, type = migration.rate_type.absolute])

   This selection policy selects n random :class:`PyGMO.individual` in the :class:`PyGMO.island`'s :class:`PyGMO.population`
   selected uniformly. This class is used exclusively in the :class:`PyGMO.island` 
   constructor as a possible kwarg for the key 's_policy'

   .. code-block:: python

      from PyGMO import *
      prob = problem.griewank(5)
      algo = algorithm.bee_colony(gen = 10) #instantiates artificial bee colony with default params and 10 generations
      random10 = migration.random_s_policy(10)
      best10 = migration.best_s_policy(10)
      isl1 = island(algo,prob,50,s_policy = best10) #10 random individuals will be selected for migration
      isl2 = island(algo,prob,50,s_policy = best10) #the 10 best individuals will be selected for migration

.. class:: PyGMO.migration.hv_greedy_s_policy([n=1, type = migration.rate_type.absolute, nadir_eps=1.0])

   This selection policy selects n :class:`PyGMO.individual` in the :class:`PyGMO.island`'s :class:`PyGMO.population`.
   The comparison between individuals is made according to the exclusive hypervolume they contribute to the population
   (see `PyGMO.hypervolume` for more details). The resulting set of individuals is created iteratively, with each step consisting of selecting the
   greatest contributor, and then removing it from the population to prevent it fromt diminishing the exclusive contributions of the
   remaining individuals.
   The reference point for all hypervolume computations is the current nadir-point of the population with an off-set determined by the ``nadir_eps``.

   NOTE: This migration applies only to multi-objective problems. In case of a single-objective problem, the `PyGMO.migration.best_s_policy` is used instead.

   .. code-block:: python

      from PyGMO import *
      prob = problem.dtlz(prob_id = 3, fdim=5)
      algo = algorithm.nsga_II(gen = 10) #instantiates the NSGA-II algorithm
      hv_greedy_10 = migration.hv_greedy_s_policy(10)
      isl = island(algo, prob, 50, s_policy = hv_greedy_10) #10 random individuals will be selected for migration

.. class:: PyGMO.migration.hv_best_s_policy([n=1, type = migration.rate_type.absolute, nadir_eps=1.0])

   This selection policy selects n :class:`PyGMO.individual` in the :class:`PyGMO.island`'s :class:`PyGMO.population`.
   The comparison between individuals is made according to the exclusive hypervolume they contribute to the population
   (see `PyGMO.hypervolume` for more details). The resulting set of individuals is created by computing all contributions
   for each of the individuals of the population, and then selecting ``n`` greatest contributors.
   The reference point for all hypervolume computations is the current nadir-point of the population with an off-set determined by the ``nadir_eps``.

   NOTE: This migration applies only to multi-objective problems. In case of a single-objective problem, the `PyGMO.migration.best_s_policy` is used instead.

   .. code-block:: python

      from PyGMO import *
      prob = problem.dtlz(prob_id = 3, fdim=5)
      algo = algorithm.nsga_II(gen = 10) #instantiates the NSGA-II algorithm
      hv_greedy_10 = migration.hv_greedy_s_policy(10)
      isl = island(algo, prob, 50, s_policy = hv_greedy_10) #10 random individuals will be selected for migration

.. class:: PyGMO.migration.fair_r_policy([n=1, type = migration.rate_type.absolute])

   A replacement policy that replaces the worst n :class:`PyGMO.individual` in the :class:`PyGMO.island`'s
   :class:`PyGMO.population` with the best n migrants. Each replacement takes place if and only if
   the migrant is considered better. If type is migration.rate_type.fractional then n, in [0,1], is interpreted
   as the fraction of the population to be replaced. This class is used exclusively in the :class:`PyGMO.island` 
   constructor as a possible kwarg for the key 'r_policy'

   .. code-block:: python

      from PyGMO import *
      prob = problem.griewank(5)
      algo = algorithm.bee_colony(10) #instantiates artificial bee colony with default params and 10 generations
      fair2 = migration.fair_r_policy(2)
      fair20pc = migration.fair_r_policy(0.2,migration.rate_type.fractional)
      isl1 = island(algo,prob,10,r_policy = fair2)  #2 of the worst individuals will be considered for replacement
      isl2 = island(algo,prob,100,r_policy = fair20pc) #20% of 100 (i.e. 20) worst individuals will be considered for replacement

.. class:: PyGMO.migration.random_r_policy([n=1, type = migration.rate_type.absolute])

   A replacement policy that replaces n random :class:`PyGMO.individual` in the :class:`PyGMO.island`'s
   :class:`PyGMO.population` with random n migrants. If type is migration.rate_type.fractional then n, in [0,1], is interpreted
   as the fraction of the population to be replaced. This class is used exclusively in the :class:`PyGMO.island` 
   constructor as a possible kwarg for the key 'r_policy'

   .. code-block:: python

      from PyGMO import *
      prob = problem.griewank(5)
      algo = algorithm.bee_colony(gen = 10) #instantiates artificial bee colony with default params and 10 generations
      random2 = migration.random_r_policy(2)
      isl = island(algo,prob,10,r_policy = random2)  #2 random individuals will be replaced with random migrants

.. class:: PyGMO.migration.worst_r_policy([n=1, type = migration.rate_type.absolute])

   A replacement policy that replaces the n worst :class:`PyGMO.individual` in the :class:`PyGMO.island`'s
   :class:`PyGMO.population` with the best n migrants. If type is migration.rate_type.fractional then n, in [0,1], is interpreted
   as the fraction of the population to be replaced. This class is used exclusively in the :class:`PyGMO.island` 
   constructor as a possible kwarg for the key 'r_policy'

   .. code-block:: python

      from PyGMO import *
      prob = problem.griewank(5)
      algo = algorithm.bee_colony(gen = 10) #instantiates artificial bee colony with default params and 10 generations
      worst2 = migration.worst_r_policy(2)
      isl = island(algo,prob,10,r_policy = worst2)  #the 2 worst individuals will be replaced by the best migrants

.. class:: PyGMO.migration.hv_greedy_r_policy([n=1, type = migration.rate_type.absolute, nadir_eps = 1.0])

   A replacement policy that replaces the worst n :class:`PyGMO.individual` in the :class:`PyGMO.island`'s
   :class:`PyGMO.population` with the best n immigrants. Each replacement takes place if and only if
   the migrant is considered better.  The comparison between individuals is made according to the exclusive hypervolume they contribute to the population
   (see `PyGMO.hypervolume` for more details).
   If type is migration.rate_type.fractional then n, in [0,1], is interpreted as the fraction of the population to be replaced.
   The reference point for all hypervolume computations is the current nadir-point of the population with an off-set determined by the ``nadir_eps``.

   The sets of worst islanders and best immigrants are determined according to the exclusive contribution by each individual.
   This is done ``iteratively``; after each request for the least (in case of worst set) or greatest (in case of best set) contributor,
   the individual is removed from the working population, in order to prevent it from diminishing the exclusive contributions by other points.

   NOTE: This migration applies only to multi-objective problems. In case of a single-objective problem, the `PyGMO.migration.fair_r_policy` is used instead.

   .. code-block:: python

      from PyGMO import *
      prob = problem.dltz3(fdim=5)
      algo = algorithm.nsga_II(gen=10) # Instantiates the NSGA-II algorithm
      hv_greedy = migration.hv_greedy_r_policy(2)
      isl = island(algo,prob,10,r_policy = hv_greedy)  # 2 of the worst individuals will be considered for replacement

.. class:: PyGMO.migration.hv_fair_r_policy([n=1, type = migration.rate_type.absolute, nadir_eps = 1.0])

   A replacement policy that replaces the worst n :class:`PyGMO.individual` in the :class:`PyGMO.island`'s
   :class:`PyGMO.population` with the best n immigrants. Each replacement takes place if and only if
   the migrant is considered better.  The comparison between individuals is made according to the exclusive hypervolume they contribute to the population
   (see `PyGMO.hypervolume` for more details).
   If type is migration.rate_type.fractional then n, in [0,1], is interpreted as the fraction of the population to be replaced.
   The reference point for all hypervolume computations is the current nadir-point of the population with an off-set determined by the ``nadir_eps``.

   The sets of worst islanders and best immigrants are determined according to the exclusive contribution by each individual.
   This is done by a single computation of the exclusive contributions over a joined population of original individuals and the immigrants.
   The vector of contributions serves as a mean for determining the required sets.

   NOTE: This migration applies only to multi-objective problems. In case of a single-objective problem, the `PyGMO.migration.fair_r_policy` is used instead.

   .. code-block:: python

      from PyGMO import *
      prob = problem.dltz3(fdim=5)
      algo = algorithm.nsga_II(gen=10) # Instantiates the NSGA-II algorithm
      hv_fair = migration.hv_fair_r_policy(3)
      isl = island(algo,prob,10,r_policy = hv_fair)  # 3 of the worst individuals will be considered for replacement

.. class:: PyGMO.distribution_type

   This class attributes are be used to set the kwarg 'distribution_type' of the :class:`PyGMO.archipelago` constructor kwarg 'migration_direction' to
   define whether the migrants will be distributed to one of the neighbouring island chosen at random or to all
   of them

   .. attribute:: PyGMO.distribution_type.point_to_point

      Migrants are distributed to one of neighbouring :class:`PyGMO.island` selected at random

   .. attribute:: PyGMO.distribution_type.broadcast

      Migrants are distributed to all neighbouring :class:`PyGMO.island` 

.. class:: PyGMO.migration_direction

   This class attributes are be used to set the kwarg 'migration_direction' of the :class:`PyGMO.archipelago` constructor kwarg 'migration_direction' to
   define whether the migrant databases will contain the incoming or the outgoing individuals.

   .. attribute:: PyGMO.migration_direction.destination

      Migrant database contains outgoing individuals

   .. attribute:: PyGMO.migration_direction.source

      Migrant database contains incoming individuals

.. class:: PyGMO.migration.rate_type

      This class attributes are used to set the second arg in the various selection and replacement policies 
      (:class:`PyGMO.migration.best_s_policy`, :class:`PyGMO.migration.fair_r_policy`, 
      :class:`PyGMO.migration.worst_r_policy`, :class:`PyGMO.migration.random_r_policy`)

   .. attribute:: PyGMO.migration.rate_type.absolute

      The number of migrants is specified as an absolute number

   .. attribute:: PyGMO.migration.rate_type.fractional

      The number of migrants is specified as fraction of the :class:`PyGMO.population` size
  
