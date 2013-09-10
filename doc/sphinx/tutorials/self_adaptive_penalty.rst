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
   In [6]: algo = algorithm.sga_gray(gen=1,cr=0.9,m=0.004,elitism=1,mutation=algorithm.sga_gray.mutation.UNIFORM,selection=algorithm.sga_gray.selection.ROULETTE,crossover=algorithm.sga_gray.crossover.SINGLE_POINT)

The algorithm has a crossover probability of 90%, a mutation probability 
of 0.4%. Furthermore, it uses a uniform mutation, a roulette selection 
and a single point crossover. For more details on this algorithm, please
refer to its documentation.

The important point here is the number of generation of the algorithm.
In fact, as the number of generation is governed by the meta-algorithm
and that the meta-algorithm modifies the fitness of the population at
each iterations, the number of generation of the algorithm must be set
to 1 to work properly.

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
(0.9999974966048448, 0.9998429119540129, 0.9999952018259526, 0.9998852014507413, 0.9999916255471595, 0.9997800290519008, 0.9999531805501336, 0.9995250701762757, 0.9999549090848836, 2.9960603414791924, 2.993714898637381, 2.9988170563822107, 0.9999217987037241)
(-14.986323158213875,)
(-0.010543942765711023, -0.0051372052770020105, -0.007791817420477187, -5.003919631359566, -5.005028396994722, -5.00114455822541, -0.0037016869694497245, -0.005798340016553993, -0.0001879930552242115)

The solution found by this algorithm is close to the optimum given with the
following:
.. code-block:: python
   In [12: print(prob.best_x)
   In [13]: print(prob.best.f)
   In [14]: print(prob.best.c)

   Out [2]:
((1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0, 3.0, 3.0, 1.0),)
((-15.0,),)
((0.0, 0.0, 0.0, -5.0, -5.0, -5.0, 0.0, 0.0, 0.0),)

Even if the solutions are really close to the optimum, the performances 
of this algorithm as described by R. Farmani and J. A. Wright in their
paper could not be retrieved with our configuration.
