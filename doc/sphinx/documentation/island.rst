Island and Archipelago
======================

NOTE: There are two different types of islands in PyGMO. The PyGMO.local_island and the PyGMO.py_island. The need for two different
types of island is purely technical (it is necessary to give the user the possibility of implementing the problem or the algorithm in python directly)
and the user needs not to know any details. We will thus here document only one class that we call
island and that, in reality, is a helper function returning automatically the correct object type.

.. class:: PyGMO.island

   At the core of the PyGMO implementation of the generalized asynchronous island model is this class: the *island*.
   An island, in PyGMO, contains a :class:`PyGMO.population` and a :class:`PyGMO.algorithm` and is the object responsible
   to launch the thread that applies the :class:`PyGMO.algorithm` to evolve the :class:`PyGMO.population`

   .. method:: PyGMO.island.__init__((PyGMO.algorithm)algo, (PyGMO.population)pop [, s_policy = best_s_policy(1), r_policy=fair_r_policy(1)])

      Constructs an island from a population. The population will be evolved by algo. Migration occurs
      at the end of each evolution if the *island* belongs to a :class:`PyGMO.archipelago`

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(2)
         algo = algorithm.de(10) #instantiates differential evolution with default params and 10 generations
         pop = population(prob,20)
         isl = island(algo,pop)

   .. method:: PyGMO.island.__init__((PyGMO.algorithm)algo, (PyGMO.problem)prob, (int)N=0 [, s_policy = best_s_policy(1), r_policy=fair_r_policy(1)])

      Constructs an island directly from a problem. The resulting population (of size N) will be evolved by algo.
      Migration occurs at the end of each evolution if the *island* belongs to an :class:`PyGMO.archipelago`

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(2)
         algo = algorithm.de(10) #instantiates differential evolution with default params and 10 generations
         isl = island(algo,prob,20)

   .. method:: PyGMO.island.evolve((int)n)

      Evolves the :class:`PyGMO.population` in the *island* performing n calls to :class:`PyGMO.algorithm`.
      At the end of each call migration occurs if the *island*  belongs to a :class:`PyGMO.archipelago`.
      Evolution happens in the background on a dedicated thread while the program flows continues.

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

   .. method:: PyGMO.island.get_evolution_time()

      Returns the time PyGMO has spent on evolving that island in milliseconds

      .. code-block:: python

         from PyGMO import *
         prob = problem.ackley(5)
         algo = algorithm.de(10) #instantiates differential evolution with default params and 10 generations
         isl = island(algo,prob,40)
         isl.evolve(100)
         isl.get_evolution_time()

   .. attribute:: PyGMO.island.algorithm

      The *island* :class:`PyGMO.algorithm`. Can be set, but not modified via its methods.

   .. attribute:: PyGMO.island.population

      The *island* :class:`PyGMO.population`. Can be set, but not modified via its methods.

   .. attribute:: PyGMO.island.problem

      A copy of the :class:`PyGMO.problem` in the :class:`PyGMO.population`. Cannot be set or modified via its methods.

.. class:: PyGMO.archipelago

   Probably the most important object in all PyGMO. An *Archipelago* is a container of :class:`PyGMO.island`
   and is responsible to start the asynchronous island model. The solutions exchange among :class:`PyGMO.island` is done
   following the routes allowed by the underlying topology

   .. method:: PyGMO.archipelago.__init__([topology = unconnected(), distribution_type = point_to_point, migration_direction = destination])

      Constructs an empty *archipelago* from a topology (defaults to :class:`PyGMO.topology.unconnected()`)
      a distribution type (defaults to :class:`PyGMO.distribution_type.point_to_point`) and a migration
      direction (defaults to :class:`PyGMO.migration_direction.destination`)

      .. code-block:: python

         from PyGMO import *
	 archi = archipelago()                            #constructs an empty archipelago with an unconnected topology
         archi = archipelago(topology = topology.ring())  #constructs an empty archipelago with a ring topology

   .. method:: PyGMO.archipelago.__init__((PyGMO.algorithm)algo, (PyGMO.problem)prob, (int)n_isl, (int)n_ind, [topology = unconnected(), distribution_type = point_to_point, migration_direction = destination])

      Constructs an empty *archipelago* from a topology (defaults to :class:`PyGMO.topology.unconnected()`)
      a distribution type (defaults to :class:`PyGMO.distribution_type.point_to_point`) and a migration
      direction (defaults to :class:`PyGMO.migration_direction.destination`).
      It then pushes back into the archipelago n_isl :class:`PyGMO.island` constructed using defaults values for the kwargs
      and (algo,prob,n_ind) as args.

      .. code-block:: python

         from PyGMO import *
	 prob = problem.grewank(30)
	 algo = algorithm.bee_colony(50)
         archi = archipelago(algo,prob,8,20)    #constructs an archipelago having 8 islands with populations
                                                #of 20 individuals. All islands have a copy of the bee_colony solver
                                                #and a copy of the griewank(30) problem

   .. method:: PyGMO.archipelago.evolve((int)n)

      Calls the method :class:`PyGMO.island.evolve` (n) on all the :class:`PyGMO.island` of the *archipelago*.
      In other words, it starts the asynchronous generalized island model that is at the core of PyGMO.

      .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(15)
         algo = algorithm.de(500)                 #instantiates differential evolution with default
                                                  #params and 100 generations
         archi = archipelago(algo,prob,8,20)
         archi.evolve(10)                         #starts the asynchronous generalized island model.
                                                  #each of the 8 islands will call algo 10 times and try to migrate in between calls

   .. method:: PyGMO.archipelago.push_back((PyGMO.island) isl)

     Pushes back isl in the archipelago taking also to also update the topological links between islands.
     This method also checks that the island inserted is compatible with the other islands present in the archipelago
     (i.e. it contains a :class`PyGMO.population` with a :class`PyGMO.problem that is compatible)

   .. method:: PyGMO.archipelago.join()

     Joins all threads (if still present) started by the  :class:`PyGMO.archipelago`. In other words it waits
     for evolution to have finished.

     .. code-block:: python

         from PyGMO import *
         prob = problem.schwefel(15)
         algo = algorithm.de(500)                 #instantiates differential evolution with default
                                                  #params and 100 generations
         archi = archipelago(algo,prob,8,20)
         archi.evolve(10)                         #starts the asynchronous generalized island model.
         archi.join()                             #waits for it to finish
	 [isl.population.champion.f for isl in archi] #builds a list with the best fittnesses found

   .. method:: PyGMO.archipelago.busy()

      Returns True if evolution is ongoing in the *archipelago*.

   .. method:: PyGMO.archipelago.interrupt()

      Halts evolution at the first occasion in all islands.

   .. automethod:: PyGMO.archipelago.draw()

   .. code-block:: python

      from PyGMO import *
      prob = problem.rosenbrock(10)
      algo = algorithm.cmaes(gen=100)
      archi = archipelago(algo,prob,8,20,topology = topology.ring())
      archi.draw()

   .. method:: PyGMO.archipelago.dump_migration_history()

      Returns a temporal history of all the archipelago migrations in one string. Each entry is in the form
      (n,src_isl,dest_isl) and logs that n individuals, from the *island* having the index src_isl,
      successfully replaced n individuals in the *population* of the *island* having the index dest_isl.


.. class:: PyGMO.zmq_island

    The ZeroMQ island behaves like any :class:`PyGMO.island` because it manages the algorithm, population and migration strategies,
    but it also communicates with other ZeroMQ islands that are in the same configuration. This is used to distribute the computation
    workload between networked machines. If it is not initialised, it will behave exactly like a local island.

   .. method:: PyGMO.zmq_island.__init__((PyGMO.algorithm)algo, (PyGMO.population)pop [, s_policy = best_s_policy(1), r_policy=fair_r_policy(1)])

        Instantiates a :class:`PyGMO.zmq_island`, with the given population and algorithm. After every evolution, the island will broadcast a copy of the population to the ZeroMQ network.

        .. code-block:: python

           from PyGMO import *
           prob = problem.schwefel(2)
           algo = algorithm.de(10)
           pop = population(prob,20)

           isl = zmq_island(algo,pop)

        Note that because the parameters specific to the ZMQ island are not set yet, this island cannot communicate and will work exactly like a local island.

   .. method:: PyGMO.zmq_island.set_broker_details((string) host, (int) port)

        Stores the broker host and port. At the moment, the broker software is a Redis server.

        .. code-block:: python

           from PyGMO import *
           isl = zmq_island(a, p)
           isl.set_broker_details("127.0.0.1", 6379)

   .. method:: PyGMO.zmq_island.set_token((string) token)

        Sets the token (also known as channel topic), which is used to determine which nodes are working on a problem with similar configurations.

        .. code-block:: python

           isl.set_token("schwefel2_de10_pop20")

        The token can be any string, but describing the algorithm and problem in the token is useful to monitor the progress of the computation. A way of monitoring the progress is described in the tutorial (TODO: link), where a user connects to the ZeroMQ network but does not perform evolution tasks.

   .. method:: PyGMO.zmq_island.set_ip((string) ip)

        This method sets the IP address to be used for incoming connections.

   .. method:: PyGMO.zmq_island.set_evolve((bool) evolve)
        
        If `set_evolve(false)` is called, then this island will not perform any evolution operations.

   .. method:: PyGMO.zmq_island.connect()

        This method will initialise the communication capabilities of the ZeroMQ island by communicating with the network, connecting to its peers and broadcasting the new connection.

        The IP address given as an argument will be used to bind the receiving socket, so it must be accessible by all of the nodes in the network. A receiving port will be chosen at random, from 1000 to 2000.

        .. code-block:: python

           from PyGMO import *
           import time

           prob = problem.schwefel(10)
           algo = algorithm.de(10)
           pop = population(prob,20)

           isl = zmq_island(algo,pop)
           isl.set_broker_details("127.0.0.1", 6379)
           isl.set_token("schwefel2_de10_pop10")
           isl.set_ip("127.0.0.1")
           isl.connect()

           while True:
               isl.evolve(10)
               print("Best: ", isl.population.champion.x)

               time.sleep(1)

        The given example will solve the 10-dimensional Schwefel problem locally, but if another instance of this program is started up in the same computer, they will exchange solutions until they converge on a solution.
