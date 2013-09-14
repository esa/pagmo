.. _immune_system:

=======================================================================
Immune system method for constraint handling
=======================================================================
The immune system method is a constraints handling method in which the 
a immune system simulation ius performed in order to reduce the number
of non feasible individuals in the population. In this tutorial, we 
will learn how to solve a constrained problem with this technique.

Method
##########
The immune system method implemented in PaGMO/PyGMO has the same 
construction as the co-evolution method. It uses two populations, 
the first population is associated to the initial problem 
from which the constraints are removed. The second population is used 
for the simulation of the immune system when antibodies are evolved 
to match a certain number of antigens. These antigens are selected to 
be the best individuals in term of feasibility in the first population.
Once found, the best antibodies are fed back into the first population.
In PaGMO/PyGMO, the matching process is done with a simple algorithm 
associated with a problem that reduce a distance between antibodies 
and antigens. This distance or matching function is simply either the
Hamming or Euclidean distance. It means that the matching process do 
not depend on the actual objective function value of the individuals. 
This method is thus computationaly efficient. The final implementation 
is based on a meta-algorithm, that takes the initial population to be 
optimized and two algorithms to evolve both the modified problem and 
the immune system. In the following we are going to see how to use 
this constraints handling technique.

Application
###########
The problem considered here is the problem g06 from the Congress on 
Evolutionary Computation 2006 (CEC2006). This problem is a cubic
function with 2 non linear inequality constraints. To solve this 
problem we will use the Differential Evolution (DE) for both the 
first and second populations. The number of iterations for the first 
algorithm must be set to 1 as the number of iterations is defined by
the meta-algorithm itself. We choose 5000 iterations here, but the
algorithm will stop before that if it has converged. The number of
iterations for the immune system is set to 70. Furthermore, we 
consider here a population size of 90 individuals.

First import the PyGMO library and choose the populations size and the
number of generation for the meta-algorithm.
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
   In [8]: pop = population(prob,pop_size);

Creates the meta-algorithm with these informations.
.. code-block:: python
   In [9]: algo_meta = algorithm.cstrs_immune_system(algorithm = algo_1, algorithm_immune = algo_2, gen = n_gen,select_method = algorithm.cstrs_immune_system.select_method.INFEASIBILITY, inject_method = algorithm.cstrs_immune_system.inject_method.CHAMPION, distance_method = algorithm.cstrs_immune_system.distance_method.EUCLIDEAN)

Here we have selected the infeasibility method where the antigen 
population is set by selecting individuals based on their 
infeasibility. The original selection where only the best individual 
is selected for the population, from COELLO did not give satisfactory 
results on this problem. The injected antibodies are a copy of the 
champion and the distance to evolve the antibodies is the Euclidean 
distance.

We can then evolve the algorithm.
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

As a comparison, you can print the best known solution for this
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
