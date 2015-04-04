Quick Start
===========

For a more advanced use of PyGMO, please refer to our :ref:`PyGMOtutorials`, 
or :ref:`PyGMOexamples`.

On one CPU
-------------------------------

Let us try to solve the 50-dimensional Schwefel problem.

.. code-block:: python

   from PyGMO import *
   prob = problem.schwefel(dim = 50)
   algo = algorithm.de(gen = 500)
   isl = island(algo,prob,20)
   print isl.population.champion.f
   isl.evolve(10)
   print isl.population.champion.f

And it is done!! We have used the algorithm Differential Evolution and we have evolved ten times 500 generations. 

.. code-block:: python

   (17643.0955597,)
   (0.0006364301698340569,)

On many CPUs ...
-----------------------------

Let us try to solve, again, the 50-dimensional Schwefel problem.

.. code-block:: python

   from PyGMO import *
   prob = problem.schwefel(dim = 50)
   algo = algorithm.de(gen = 500)
   archi = archipelago(algo,prob,8,20)
   print min([isl.population.champion.f for isl in archi])
   archi.evolve(10)
   print min([isl.population.champion.f for isl in archi])

And it is done!! We have launched eight separated threads each one running an instance of Differential Evolution.
Each thread evolves separately ten times 500 generations. We then print the best found in the 8 runs. 

... and migrating solutions ...
----------------------------------------------------

Let us try to solve, again, the 50-dimensional Schwefel problem.

.. code-block:: python

   from PyGMO import *
   prob = problem.schwefel(dim = 50)
   algo = algorithm.de(gen = 500)
   archi = archipelago(algo,prob,8,20, topo = topology.ring())
   print min([isl.population.champion.f for isl in archi])
   archi.evolve(10)
   print min([isl.population.champion.f for isl in archi])

And it is done!! We have launched eight separated threads each one running an instance of Differential Evolution. Each thread evolves
for 500 generation its population, then it exchanges solutions according to the defined ring topology. All happens asynchronously in the background
We then print the best found in the 8 runs. 


... and between different algorithms ...
-------------------------------------------------------------------------

Let us try to solve, again, the 50-dimensional Schwefel problem.

.. code-block:: python

   from PyGMO import *
   prob = problem.schwefel(dim = 50)
   algo = []
   for i in range(1,9):
       algo.append(algorithm.de(gen=500,variant=i))
   archi = archipelago(topo=topology.ring())
   for i in range(0,8):
       archi.push_back(island(algo[i],prob,20))
   print min([isl.population.champion.f for isl in archi])
   archi.evolve(20)
   print min([isl.population.champion.f for isl in archi])

And it is done!! We have instantiated 8 different variants of Differential Evolution
and have them cooperate to solve the same optimization problem!

