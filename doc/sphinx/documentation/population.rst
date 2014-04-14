Population
============

.. class:: PyGMO.population

   This class represents the concept of a *population* as a collection of :class:`PyGMO.individual` evaluated w.r.t. a :class:`PyGMO.problem`.
   A *population* in PyGMO makes sure that all the individuals it contains are consistent w.r.t. the problem and it constantly keeps
   its :class:`PyGMO.champion` updated. Also (relevant only for multi-objective problems) the *population*, keeps constantly updated
   1) a domination list that is a list containing, per individual I, the indexes of the individuals I dominates, and 2) a domination count
   that is a list containing, per individual I, the number of individuals that dominate individual I.
   From an evolutionary point of view one can see the *population* as a set of individuals living in an environment (the :class:`problem`) which defines their fitness values

   .. method:: __init__((PyGMO.problem)prob [, (int)n_individuals, (int)seed])

      Constructs a population containing n_individuals (default is an empty population) evaluated w.r.t. prob. 
      Each individual gets initialized at random with his chromosome in [:class:`problem.lb`, :class:`problem.ub`] and his velocity 
      in [(:class:`problem.lb`-:class:`problem.ub`)/2, (:class:`problem.ub`- :class:`problem.lb`)/2].
      The initialization process of the population can be seeded by providing the argument `seed`.

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(50)
         pop1 = population(prob) #empty population
         pop2 = population(prob,20) #population with 20 individuals
         pop3 = population(prob,20,123) #population with 20 individuals intialized with seed=123

   .. method:: __init__((PyGMO.population)pop)

      Constructs a new *population* performing a deep-copy of a second population pop. 

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(50)
         pop1 = population(prob,10) #population with 10 individuals
         pop2 = population(pop1)    #a second identical population with 10 individuals 
                                    #(the problem and the individuals are copies of the original pop1)

   .. method:: push_back((list) x)

      Appends an :class:`individual` having chromosme x to the population, if compatible with the :class:`problem`. Its velocity
      is initialized at random, its memory is set equal to the current position.

      NOTE: There is no way to push_back into a *population* directly an :class:`PyGMO.individual`
      This design choice was made as there would be no way to ensure efficiently that the information contained
      in the :class:`PyGMO.individual` is consistent (if not by wasting one function evaluation) with the :class:`PyGMO.problem` defining the *population*

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(2)
         pop = population(prob)
         pop.push_back([1.12,2.34])

   .. method:: erase((int) idx)

      Erases the :class:`individual` with index idx from the *population*. Domination list and 
      domination count are updated accordingly.

      NOTE: after such an operation all indexes will be renamed so that if the individual with idx = n is erased, 
      after the erase has completed the individual that had idx=n+1 will have idx = n

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(2)
         pop = population(prob,5)
         pop.erase(0)

   .. method:: set_x((int)idx, (list) x)

      Sets the chromosome of the :class:`PyGMO.individual` with index idx in the population to x. Updates autatically the memory and 
      the *population* :class:`champion`

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(2)
         pop = population(prob,2)
         pop.set_x(0,[3.12,4.56])

   .. method:: set_v((int)idx, (list) v)

      Sets the velocity of the :class:`PyGMO.individual` with index idx in the population to v

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(2)
         pop = population(prob,2)
         pop.set_v(0,[0.12,-0.22])

   .. method:: get_domination_list((int)idx)

      Returns a list containing all the indexes of the individual dominated by the individual with index idx

      .. code-block:: python

         from PyGMO import *
         prob = problem.zdt(1)
         pop = population(prob,10)
         ls = pop.get_domination_list(1)

   .. method:: get_domination_count((int)idx)

      Returns the domination count for the individual idx (that is how many individuals in the population dominate idx?)

      .. code-block:: python

         from PyGMO import *
         prob = problem.zdt(1)
         pop = population(prob,10)
         c = pop.get_domination_count(1)

   .. method:: compute_pareto_fronts()

      Returns the Pareto fronts of the population in form of a list of lists each one containing the idx
      of the individuals belonging to a particular Pareto Front

      .. code-block:: python

         from PyGMO import *
         prob = problem.zdt(1)
         pop = population(prob,10)
         pf = pop.compute_pareto_fronts()

   .. method:: plot_pareto_fronts(comp = [0,1])

      Plots the pareto fronts in a sliced 2-D graph representing the two objective function components specified
      in comp

      .. code-block:: python

         from PyGMO import *
         prob = problem.zdt(1)
         pop = population(prob,100)
         pf = pop.plot_pareto_fronts()

   .. method:: get_best_idx((int) n)

      Returns the n best indexes of the :class:`PyGMO.individual` in a *population*. The best 
      :class:`PyGMO.individual` are computed according to non-dominated sorting in populations that
      have a multi-objective problem.

   .. method:: get_worst_idx()

      Returns the index of the worst :class:`PyGMO.individual` in a *population*. The worst 
      :class:`PyGMO.individual` is computed according to non-dominated sorting in populations that
      have a multi-objective problem.

      .. code-block:: python

         from PyGMO import *
         prob = problem.zdt(3)
         pop = population(prob,3) #population with 3 individuals
         best_guy = pop.get_best_idx()
         worst_guy = pop.get_worst_idx()

   .. method:: mean_velocity()
 
      Evaluates the *population* mean velocity

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(3)
         pop = population(prob,30) 
         v = pop.mean_velocity()

   .. method:: race((int) n_winners, (int) min_trials=0, (int) max_feval=500, (float) delta=0.05, (list) racers_idx=[])

	  Races individuals in a population

	  * n_winners: number of winners in the race
	  * min_trials: minimum amount of evaluations before an individual can stop racing
	  * delta: Statistical test confidence
	  * racers_idx: indices of the individuals in pop to be raced

   .. method:: repair((int) idx, (:class:`problem`) repair_algo)

	  Repairs the individual at the position idx

	  * idx: index of the individual to repair
	  * repair_algo: 'repairing' optimization algorithm to use. It should be able to deal with population of size 1.

   .. attribute:: champion
      :noindex:

      Returns a copy of the *population* :class:`PyGMO.champion` 

   .. attribute:: problem
      :noindex:

      Returns a copy of the :class:`problem` in the *population*

      NOTE: since it is only a copy that is returned, it is impossible to modify a problem in a population
      directly. The following code is thus WRONG as it changes the bounds of an instance of the problem that is created
      on the fly and then destroyed

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(3)
         pop = population(prob,30) 
         lb = list(prob.lb)
         ub = list(prob.ub)
         lb[0]=-10
         pop.problem.set_bounds(lb,ub) #This line is completely uneffective ...


