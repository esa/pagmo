.. _self_adaptive_penalty:

=======================================================================
Self-adaptive penalty
=======================================================================
The self-adaptive penalty is a method to handle constraints in
optimization problems. In this tutorial, we are going to solve a
constrained problem using the self-adaptive technique.

Method
##########
The self-adaptive is a part of what is called penalty constrained
handling techniques. The idea is to penalize in term of cost function
the individuals that do not satisfy the constraints. Unfortunately,
the penalization must be well tuned for the method to be effective.
Given a population, the self-adaptive technique automatically adapt
the penalty to apply, depending on some infeasibility measure of
the individuals. The aim is not to reject completely indifiduals that 
are not totally feasible, or in other words that slightly violate some
of the constraints. The actual implementation of the self-adaptive 
method in PaGMO/PyGMO is based on population fitness evaluation 
modification.

Application
###########
We consider here the constrained problem g01 from the Congress on 
Evolutionary Computation 2006 (CEC2006). This problem is a quadratic
problem with 9 linear inequality constraints and 6 active constraints.
To solve it, we use a simple genetic algorithm with gray encoding
coupled with this constraints handling technique. To do so in 
PaGMO/PyGMO, we use a meta-algorithm that modifies the original 
algorithm that was only dedicated to solve unconstrained problems. In 
the population, 70 individuals are present. 

Step by step, we first import the PyGMO library and choose the
populations size and the number of generation for the meta-algorithm.

.. code-block:: python

   In [1]: from PyGMO import *
   In [2]: pop_size = 70
   In [3]: n_gen = 20000

Then we create the problem, the associated population and the algorithm.

.. code-block:: python

   In [4]: prob = problem.cec2006(1)
   In [5]: pop = population(prob,pop_size)
   In [6]: algo = algorithm.sga_gray(gen=1,cr=0.9,m=0.004,elitism=0,mutation=algorithm.sga_gray.mutation.UNIFORM,selection=algorithm.sga_gray.selection.ROULETTE,crossover=algorithm.sga_gray.crossover.SINGLE_POINT)

The algorithm has a crossover probability of 90%, a mutation probability 
of 0.4%. Furthermore, it uses a uniform mutation, a roulette selection 
and a single point crossover. For more details on this algorithm, please
refer to its documentation.

The important point here is the number of generation of the algorithm.
In fact, as the number of generation is governed by the meta-algorithm
and that the meta-algorithm modifies the fitness of the population at
each iterations, the number of generation of the algorithm must be set
to 1 to work properly. The same comment holds for the elitism, with
this algorithm, it is possible to 

The next step consists in creating the meta-algorithm.

.. code-block:: python

   In [7]: algo_self_adaptive = algorithm.cstrs_self_adaptive(algo, n_gen)

Then we evolve it and print the solutions.

.. code-block:: python

   In [8]: pop = algo_self_adaptive.evolve(pop)
   In [9]: print(pop.champion.x)
   In [10]: print(pop.champion.f)
   In [11]: print(pop.champion.c)
   Out [1]:
   (0.9999999701976767, 0.9999997615814138, 0.9999998509883836, 1.0, 1.0, 0.9999980032443405, 0.9999997615814138, 0.9999997615814138, 1.0, 2.9999972283839353, 2.99999126791928, 2.99997934698997, 0.9999922513959483)
   (-14.999955534934077,)
   (-0.008994877606477658, -0.0184053188087141, -0.01848650033731758, -5.004269987471997, -5.004496604338187, -4.982709675512006, -0.0023840070481302433, -0.004443645609725877, -0.0037382544201092216)

The solution found by this algorithm is close to the optimum given
with the following:

.. code-block:: python

   In [12: print(prob.best_x)
   In [13]: print(prob.best.f)
   In [14]: print(prob.best.c)
   Out [2]:
   ((1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0, 3.0, 3.0, 1.0),)
   ((-15.0,),)
   ((0.0, 0.0, 0.0, -5.0, -5.0, -5.0, 0.0, 0.0, 0.0),)

Even if the solutions are really close to the optimum, the exact 
same performances of this algorithm as described by R. Farmani 
and J. A. Wright in their paper could not be retrieved with our 
configuration.
