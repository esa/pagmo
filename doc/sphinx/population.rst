A Population
============

.. class:: population

   This class represents the concept of a *population* as a collection of :class:`individual` evaluated w.r.t. a :class:`problem`.
   A *population* in PyGMO makes sure that all the individuals it contains are consistent w.r.t. the problem and it constantly keeps
   its :class:`champion` updated. Also, the *population*, keeps constantly updated the domination list, that
   is a list containing, per individual I, the individuals that I dominates (only useful for multiobjective optimization)
   From an evolutionary point of view one can see the *population* as a set of individuals
   living in an environment (the :class:`problem`) which defines their fitness values

   .. method:: __init__((PyGMO.problem)prob [, (int)n_individuals])

      Constructs a population containing n_individuals (default is an empty population) evaluated w.r.t. prob. 
      Each individual gets initialized at random with his chromosome in [:class:`problem.lb`, :class:`problem.ub`] and his velocity 
      in [(:class:`problem.lb`-:class:`problem.ub`)/2, (:class:`problem.ub`- :class:`problem.lb`)/2]

      .. code-block:: python


         from PyGMO import *
         prob = problem.schwefel(50)
         pop1 = population(prob) #empty population
         pop2 = population(prob,20) #population with 20 individuals

   .. method:: __init__((PyGMO.population)pop)

      Constructs a new *population* performing a deep-copy of a second population pop. 

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(50)
         pop1 = population(prob,10) #population with 10 individuals
         pop2 = population(pop1)    #a second identical population with 10 individuals 
                                    #(the problem and the individuals are copies of the original pop1)

   .. method:: push_back((list) x)

      Appends the :class:`individual` having chromosme x to the population, if compatible with the :class:`problem`. Its velocity
      is initialized at random, its memory is set equal to the current position.

      NOTE: There is no way to push_back into a *population* directly an :class:`individual`.
      This design choice was made as there would be no way to ensure efficiently that the information contained
      in the :class:`individual` is consistent (if not by wasting one function evaluation) with the :class:`problem` defining the *population*

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(2)
         pop = population(prob)
         pop.push_back([1.12,2.34])

   .. method:: set_x((int)idx, (list) x)

      Sets the chromosome of the :class:`individual` with index idx in the population to x. Updates autatically the memory and 
      the *population* :class:`champion`

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(2)
         pop = population(prob,2)
         pop.set_x(0,[3.12,4.56])

   .. method:: set_v((int)idx, (list) v)

      Sets the velocity of the :class:`individual` with index idx in the population to v

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(2)
         pop = population(prob,2)
         pop.set_v(0,[0.12,-0.22])

   .. method:: get_domination_list((int)idx)

      Returns a list containing all the indexes of the individual dominated by the individual with index idx

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(2)
         pop = population(prob,10)
         pop.get_domination_list(1)

   .. method:: get_best_idx()

      Returns the index of the best :class:`individual` in a *population*. The best :class:`individual` is the one dominating the most
      number of individuals

   .. method:: get_worst_idx()

      Returns the index of the worst :class:`individual` in a *population*. The best :class:`individual` is the one dominating the least
      number of individuals

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(3)
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

   .. attribute:: champion

      Returns a copy of the *population* :class:`champion` 

   .. attribute:: problem

      Returns a copy of the :class:`problem` in the *population*

      NOTE: since it is only a copy that is returned it is impossible to set the bounds of a problem in a population
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


 


