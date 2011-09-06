Island and Archipelago
======================

NOTE: There are two different types of islands in PyGMO. The PyGMO.local_island and the PyGMO.py_island. The need for two different
types of island is purely technical (it is necessary to give the user the possibility of implementing the problem or the algorithm in python directly)
and the user needs not to know any details. We will thus here document only one class that we call
island and that, in reality, is a helper function returning automatically the correct object type. 

.. class:: PyGMO.island()

   At the core of the PyGMO implementation of the generalized asynchronous island model is this class: the *island*.
   An island, in PyGMO, contains a :class:`PyGMO.population` and a :class:`PyGMO.algorithm` and is the object responsible
   to launch the thread that applies the :class:`PyGMO.algorithm` to evolve the :class:`PyGMO.population`

   .. method:: PyGMO.island.__init__((PyGMO.algorithm)algo, (PyGMO.population)pop[ , migr_prob=1., s_policy = PyGMO.migration., r_policy=])

      Constructs an island from a population. The population will be evolved by algo. Migration occurs with
      probability migr_prob at the end of each evolution if the *island* belongs to a :class:`PyGMO.archipelago`

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(2)
         algo = algorithm.de(10) #instantiates differential evolution with default params and 10 generations
         pop = population(prob,20)
         isl = island(algo,pop)

   .. method:: PyGMO.island.__init__((PyGMO.algorithm)algo, (PyGMO.problem)prob, (int)N[ , migr_prob=1.])

      Constructs an island directly from a problem. The resulting population (of size N) will be evolved by algo. 
      Migration occurs with probability migr_prob at the end of each evolution if the *island* belongs to an :class:`PyGMO.archipelago`

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(2)
         algo = algorithm.de(10) #instantiates differential evolution with default params and 10 generations
         isl = island(algo,prob,20)

   .. method:: PyGMO.island.evolve(int(n))
  
      Evolves the :class:`PyGMO.population` in the *island* performing n calls to :class:`PyGMO.algorithm`.
      At the end of each call migration occur with probability migr_prob if the *island*  belongs to a :class:`PyGMO.archipelago`.
      Evolution happens in the background on a new thread while the program flows continues.
      
      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(1500)
         algo = algorithm.de(100)       #instantiates differential evolution with default params and 100 generations
         pop = population(prob,20)  
         isl = island(algo,pop)
         isl.evolve(1)                  #calls algo to evolve pop. The program flow continues

   .. method:: PyGMO.island.join()

      Joins the thread (if still present) started by the  :class:`PyGMO.island.evolve` method. In other words it waits
      for :class:`PyGMO.island.evolve` to have finished.

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(1500)
         algo = algorithm.de(100) #instantiates differential evolution with default params and 100 generations
         pop = population(prob,20)  
         isl = island(algo,pop)
         isl.evolve(1) #calls algo to evolve pop. The program flow continues
         isl.join() #Waits for the evolve to finish (i.e. waits for completion of the 100 generations of differential evolution)

   .. method:: PyGMO.island.busy()

         Returns True if evolution is ongoing in the *island*.

   .. method:: PyGMO.island.set_x((int)idx,(list) x)

      Sets a new chromosome for the idx-th :class:`PyGMO.individual` in the :class:`PyGMO.population` 
      of the *island* to x.

      .. code-block:: python

         from PyGMO import *
         prob = problem.ackley(5)
         algo = algorithm.de(10)             #instantiates differential evolution with default params and 10 generations
         isl = island(algo,prob,10)
	 isl.population.set_x(0,[1,2,3,4,5]) # This is completely uneffective 
                                             # as the 'attribute' population returns a copy
         isl.set_x(0,[1,2,3,4,5])            # This works!!
  

   .. method:: PyGMO.island.set_v((int)idx,(list) v)

      Sets the velocity of the idx-th :class:`PyGMO.individual` in the :class:`PyGMO.population` 
      of the *island* to v.

      .. code-block:: python

         from PyGMO import *
         prob = problem.ackley(5)
         algo = algorithm.de(10) #instantiates differential evolution with default params and 10 generations
         isl = island(algo,prob,10)
	 isl.population.set_v(0,[0.02,0.03,-0.3,0.12,0.1]) # This is completely uneffective 
                                                           # as the 'attribute' population returns a copy
         isl.set_v(0,[0.02,0.03,-0.3,0.12,0.1])            # This works!!

   .. attribute:: PyGMO.island.algorithm

      The *island* :class:`PyGMO.algorithm`. Can be set, but not modified via its methods.

   .. attribute:: PyGMO.island.population

      The *island* :class:`PyGMO.population`. Can be set, but not modified via its methods.

   .. attribute:: PyGMO.island.problem

      A copy of the :class:`PyGMO.problem` in the :class:`PyGMO.population`. Cannot be set or modified via its methods.

   .. attribute:: PyGMO.island.migration_probability

      The probability that migration occur at the end of each evolution. Cannot be set (maybe in future versions)



