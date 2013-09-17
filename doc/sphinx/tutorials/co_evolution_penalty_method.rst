.. _co_evolution_penalty_method:

=======================================================================
Co-evolution penalty method
=======================================================================

Co-evolution mimics biology where distinct populations depends on each
other for their evolution. This technique has been used by some authors
to deal with constraints in optimization and has been implemented in
PaGMO/PyGMO. In this tutorial, we will learn how to solve a constrained
problem with this technique.

Method
##########
The co-evolution method implemented in PaGMO/PyGMO has two populations
evolving at the same time. The first one is a population associated 
with the initial problem from which the constraints are removed and
where penalties are applied for individuals that do not satisfy the
constraints. To apply the penalty a certain weight is needed to
penalize more or less the individuals in order to drive the
optimizer towards the feasible domain. As the weights are cumbersome to
determine, the co-evolution updates the weights by the mean of a second
population for which the fitness depends on the constraints
satisfaction of the first population. In PaGMO/PyGMO this method is
implemented within a meta-algorithm. The meta-algorithm takes the 
initial population to be optimized and two algorithms to evolve
both the first and second populations. In the following we are going
to see how to effectively solve a constrained problem.

Application
###########
The problem considered here is the problem g05 from the Congress on 
Evolutionary Computation 2006 (CEC2006). This problem is a cubic
function with 2 linear inequality and 3 non linear equality 
constraints. To solve this problem we will use the Differential
Evolution (DE) for both the first and second populations. The 
first population size is set to 60 and the second to 30. In practice
it is recommended to take the second population size as being smaller
than the first population. Half the size of it seems to be a good
compromize.

For the choice of the number of generations, we are going to evolve
the first algorithm for 25 generations. The second algorithm number
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

   In [5]: algo_1 = algorithm.de(gen = n_gen, xtol=1e-30, ftol=1e-30)
   In [6]: algo_2 = algorithm.de(gen = 1, xtol=1e-30, ftol=1e-30)

Select the problem and associate a population to this problem.

.. code-block:: python

   In [7]: prob = problem.cec2006(5)
   In [8]: pop = population(prob,pop_1_size);

Creates the meta-algorithm with these informations.

.. code-block:: python

   In [9]: algo_coevo = algorithm.cstrs_co_evolution(algorithm = algo_1, algorithm_2 = algo_2, pop_2_size = pop_2_size, gen = n_gen)

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

The solution found by this method is the optimum of the constrainted
problem. In practice, to get this result, you might need to run
this tutorial more than once as with the differential evolution
algorithm, the meta-algorithm can be trapped in a local optimum.
