.. _co_evolution_penalty_method:

=======================================================================
Co-evolution penalty method
=======================================================================

Co-evolution mimics the biological behaviour of two distinct populations 
whose evolutions are interdependent. This technique has been used by some authors
to deal with constraints in optimization and has been implemented in
PaGMO/PyGMO. In this tutorial, we will learn how to solve a constrained
problem with this technique.

Method
##########
The co-evolution method implemented in PaGMO/PyGMO has two populations
evolving at the same time. The first one is a population associated 
to the initial problem from which the constraints are removed and
where a penalized fitness is assigned to individuals that do not satisfy the
constraints. To apply such a penalization the penalty coefficients need to
be assessed, according to the level of infeasibility of the individuals in the population, 
in order to drive it towards feasible region of the search space. 
As the coefficients are cumbersome to determine, the co-evolution technique 
updates them by mean of a second population that encodes them. 
The fitness of the second population depends on the constraints
satisfaction of the individuals in the first population. 
Hence the co-evolution happens sequentially: a set of penalty coefficents are firstly 
assigned ramdomly (the individuals in the second population). For each set of coefficient 
(for each individual in the second population) the first population is evolved 
for a certain number of iterations. Then the second population is evolved, with a 
fitness definition that depends on the results of the evolution of the first population 
and process then iterates with the new set of penalty coefficients.

In PaGMO/PyGMO this method is implemented within a meta-algorithm. The meta-algorithm takes the 
initial population to be optimized and two algorithms to evolve
respectively the first and second populations. In the following we are going
to see how to effectively solve a constrained problem.

Application
###########
The problem considered here is the problem g05 from the Congress on 
Evolutionary Computation 2006 (CEC2006). This problem has a cubic
objective function with two linear inequality constraints and 3 non linear equality 
constraints. 
To solve this problem the Differential
Evolution (DE) is chosen for both the first and second population. The 
first population size is set to 60 and the second to 30. It is good practice to set
the second population size smaller than the first population size. 
Half the size of it seems to be experimentally a good compromize.
The first algorithm evolves for 25 generations while the second algorithm number
of generations must be set to 1. This is important as the second
population uses the information from the first population and
can not be evolved more than one time. This is repeated for 50 
generations until convergence occurs.

First import the PyGMO library and choose the populations size and the
number of generation for the meta-algorithm.

.. code-block:: python
   
   In [1]: from PyGMO import *
   In [2]: pop_1_size = 60
   In [3]: pop_2_size = 30
   In [4]: n_gen = 50

Then creates the algorithms you wish to use for both populations. Here
we have decided to use the Differential Evolution for both populations.

.. code-block:: python

   In [5]: algo_1 = algorithm.de(gen = 25, xtol=1e-30, ftol=1e-30)
   In [6]: algo_2 = algorithm.de(gen = 1, xtol=1e-30, ftol=1e-30)

Select the problem and associate a population to this problem.

.. code-block:: python

   In [7]: prob = problem.cec2006(5)
   In [8]: pop = population(prob,pop_1_size)

Creates the meta-algorithm with these informations.

.. code-block:: python

   In [9]: algo_coevo = algorithm.cstrs_co_evolution(original_algo = algo_1, original_algo_penalties = algo_2, pop_penalties_size = pop_2_size, gen = n_gen)

Evolve the algorithm.

.. code-block:: python

   In [10]: pop = algo_coevo.evolve(pop)

And finally, print the solutions.

.. code-block:: python

   In [11]: print(pop.champion.x)
   In [12]: print(pop.champion.f)
   In [13]: print(pop.champion.c)

   Out [1]:
   (679.9451523687442, 1026.0669716493055, 0.11887636619084481, -0.3962334865933514)
   (5126.4967140070985,)
   (9.999999997489795e-05, 9.999999997489795e-05, 9.999999997489795e-05, -0.03489014721580386, -1.0651098527841962)

The solution found by this method is the global optimum of the constrainted
problem. Due to the stochastic behavior of the algorithm performing multiple 
runs is always raccomanded.
