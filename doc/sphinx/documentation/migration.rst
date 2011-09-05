Migration 
=========

The algorithm
-------------

Migration, in PyGMO, happens asynchronously in each :class:`PyGMO.island` between calls of the evolve() method 
and only when the islands are insterted in a :class:`PyGMO.archipelago`. The 
algorithm is rather complex and the user does not need to know/understand its details as PyGMO sets defaults values for all
of its many parameters. These are:

* Migration Rate - Default is 1

 This parameter defines the number of individuals that are selected from each island for migration. 
 To Migration Rate is set by the 's_policy' and 'r_policy' kwarg in the island constructor. This can be done by specifying
 the number of individuals that need to migrate (migration.rate_type.absolute) or the fraction of the population individuals
 (migration.rate_type.fractional)

 .. code-block:: python

      from PyGMO import *
      prob = problem.schwefel(15)
      algo = algorithm.de(100) #instantiates differential evolution with default params and 100 generations
      selection = migration.best_s_policy(0.25,migration.rate_type.fractional)
      replacement = migration.fair_r_policy(0.25,migration.rate_type.fractional)
      isl = island(algo,prob,s_policy = selection, r_policy = replacement)

* Migration Direction - Default is 'destination'

 In PyGMO the asynchronous migration is implemented by keeping a migration database on each island. Then one of the following
 options can be followed:

 * 'destination': The internal migration database stores, for each island, the individuals that are meant to migrate
   from that island. Before each evolution, the island will get migrating individuals from those made available
   by the islands connecting to it. After each evolution, the island will update its list of best individuals in the database.

 * 'source': The internal migration database stores for each island the individuals that are meant to migrate
   to that island. Before each evolution, an island will check if individuals destined to it are available in the database,
   and, in such case will, migrate over incoming individuals before starting evolution.
   After each evolution, the island will place its candidate individuals for emigration in the database slots of the island(s) to which
   it connects.

 The migration direction is set by the 'migration_direction' kwarg in the archipelago constructor

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(15)
         algo = algorithm.de(100) #instantiates differential evolution with default params and 100 generations
         direction = migration_direction.source
         archi = archipelago(migration_direction = direction)
 
* Migration Distribution Type - Default is 'point to point'

 When migration happens one has to decide which of the connected islands contributes to the event. This is decided
 by the distribution type that can be one of the following:

 * 'point to point': only one of the neighbourghing islands, selected at random, is sending (or receiving) the individuals

 * 'broadcast': all neighbourghing islands are sending (or receiving) the individuals

 The migration distribution type is set by the 'distribution_type' kwarg in the archipelago constructor

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(15)
         algo = algorithm.de(100) #instantiates differential evolution with default params and 100 generations
         distribution = distribution_type.broadcast
         archi = archipelago(distribution_type = distribution)

* Migration Selection Policy - Default is 'best_s_policy'

 The selection policy is the object responsible to choose out of a population the individuals that will migrate. All
 selection policies derive from the same base class and currently a few are implemented:

 * 'migration.best_s_policy': simply selects the best individuals

 The selection policy is set by the 's_policy' kwarg in the island constructor

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(15)
         algo = algorithm.de(100) #instantiates differential evolution with default params and 100 generations
         best2 = migration.best_s_policy(2) # two individuals will be selected as the best
         isl = island(algo,prob,s_policy = best2)

* Migration Replacement Policy - Default is 'fair_r_policy'

 The replacement policy is the object responsible to substitute the individuals in a population with the
 migrants. All replacement policies derive from the same base class and currently a few are implemented:

 * 'migration.fair_r_policy': simply replaces the worst individuals in the island  with the best of the incoming migrants. This is subject to the added condition that the migrants are better.

 * 'migration.random_r_policy': replaces random individuals in the island with random incoming migrants

 * 'migration.worst_r_policy': replaces the worst individuals in the island with the best of the incoming migrants.

 The replacement policy is set by the 'r_policy' kwarg in the island constructor

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(15)
         algo = algorithm.de(100) #instantiates differential evolution with default params and 100 generations
         random4 = migration.random_policy(4) # four individuals will be selected at random 
					      # from the migrants and will replace random individuals
         isl = island(algo,prob,s_policy = best2)


* Migration Probability - Default is 1

 The migration probability determines whether migration occurs at all between calls of the evolve() method. 
 It is set by the 'migr_prob' kwarg of the island constructor.

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(15)
         algo = algorithm.de(100) #instantiates differential evolution with default params and 100 generations
         p = 0.5
         isl = island(algo,prob,migr_prob = p)

* Migration Topology (i.e. which island is connected to which island) - Default is 'unconnected'


The Classes
---------------------

  