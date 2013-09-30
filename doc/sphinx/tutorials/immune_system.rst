.. _immune_system:

=======================================================================
Immune system method for constraint handling
=======================================================================
The immune system method is a constraints handling method in which  
a immune system simulation is performed in order to reduce the number
of non feasible individuals in the population through the action of an 
antibody population. In this tutorial, we 
will learn how to solve a constrained problem with this technique.

Method
##########
The immune system method implemented in PaGMO/PyGMO has a similar 
design to the one of the co-evolution method. It uses two populations, 
the first population is associated to the initial problem 
from which the constraints are removed. The second population emulates 
the immune system where antibodies are evolved 
to match a certain number of antigens. These antigens are selected to 
be the best individuals, in term of feasibility, of the first population.
Once found, the best antibodies are fed back into the first population.
In PaGMO/PyGMO, the matching process is done with a simple algorithm 
associated with a problem that reduces the distance between antibodies 
and antigens. This distance or matching function is either the
Hamming or Euclidean distance. It means that the matching process does 
not require additional evaluations of the objective or constraints functions. 
This step is thus computationaly efficient. The final implementation 
is based on a meta-algorithm, that takes the initial population to be 
optimized and two algorithms to evolve both the population associated to
the modified problem and the immune system. 
In the following we are going to see how to use this constraints handling technique.

Application
###########
The problem considered here is the problem g06 from the Congress on 
Evolutionary Computation 2006 (CEC2006). 
This problem has a cubic
objective function with two non linear inequality constraints. 
To solve this problem the Differential Evolution (DE) algorithm is used for both the 
first and the second population. The number of iterations for the first 
algorithm must be set to 1 as the overall number of iterations is driven by
the meta-algorithm itself. The number of iteration of the meta-algorithm is set to 5000, 
however the algorithm will stop as soon as it reaches convergence. The number of
iterations for the immune system is set to 70 and a initial population of 90 individuals is chosen.

First import the PyGMO library and choose the population size and the
number of generations for the meta-algorithm.

.. code-block:: python
   
   In [1]: from PyGMO import *
   In [2]: pop_size = 90
   In [3]: n_gen = 5000
   In [4]: n_immune_gen = 70

Then creates the algorithms you wish to use for both populations. Here
we have decided to use the Differential Evolution for both populations.

.. code-block:: python

   In [5]: algo_1 = algorithm.de(gen = 1, xtol=1e-30, ftol=1e-30)
   In [6]: algo_2 = algorithm.de(gen = n_immune_gen, xtol=1e-30, ftol=1e-30)

Select the problem and associate a population to this problem.

.. code-block:: python

   In [7]: prob = problem.cec2006(6)
   In [8]: pop = population(prob,pop_size)

Creates the meta-algorithm with these informations.

.. code-block:: python

   In [9]: algo_meta = algorithm.cstrs_immune_system(algorithm = algo_1, algorithm_immune = algo_2, gen = n_gen,select_method = algorithm.cstrs_immune_system.select_method.INFEASIBILITY, inject_method = algorithm.cstrs_immune_system.inject_method.CHAMPION, distance_method = algorithm.cstrs_immune_system.distance_method.EUCLIDEAN)

Here we have selected the infeasibility method where the antigen 
population is set by selecting individuals based on their 
infeasibility. The original selection, where only the best infeasible individual 
is selected for the population, from Coello did not give satisfactory 
results on this problem. The injected antibodies are a copy of the 
champion and the distance to evolve the antibodies is the Euclidean 
distance.

Evolve the population with the mat-algorithm described.

.. code-block:: python

   In [10]: pop = algo_meta.evolve(pop)

And finally, print the solutions.

.. code-block:: python

   In [11]: print(pop.champion.x)
   In [12]: print(pop.champion.f)
   In [13]: print(pop.champion.c)
   Out [1]:
   (14.094999999999994, 0.8429607892154646)
   (-6961.813875580156,)
   (-3.552713678800501e-15, 0.0)

As a comparison, the best known solution can be printed for this
particular problem:

.. code-block:: python

   In [11]: print(prob.best_x)
   In [12]: print(prob.best_f)
   In [13]: print(prob.best_c)
   Out [2]:
   ((14.095, 0.8429607892154796),)
   ((-6961.813875580138,),)
   ((-7.105427357601002e-15, 0.0),)

As seen, the algorithm has converged to the optimal constrained 
solution.
