.. _self_adaptive_penalty:

=======================================================================
Self-adaptive penalty
=======================================================================
The self-adaptive penalty is a method to handle constraints in
optimization problems. It is an adaptive penalty method where the 
penalization factor is computed dynamically according to the 
individuals in the population at the current generation.
In this tutorial, it is explained how to solve a single objective
constrained problem using the self-adaptive technique.

Method
##########
The self-adaptive method is a part of what is called adaptive penalty
methods for constrained optimization. 
The idea is to penalize in term of cost function
the infeasible individuals. Minimization is assumed in the design of the algorithm.

The infeasibility measure is computed through a two-stage penalty approach. 
The first stage applies just if there is one or more infeasible solution in the 
population that has a lower, hence potentially better, objective function value 
than the best feasible solution in the population. The penalty is hence applied 
to all the infeasible solutions in the population according to the relative 
order of infeasibility, in such a way that the solution with a high infeasibility 
rate but a low objective function has a fitness value close to the one of the best 
individual in the population. This is to avoid the elimination of infeasible solutions 
potentially close to the real optimum.
The second stage penalizes exponentially all the infeasible individual according to 
their infeasibility rate and a scaling factor that is a measure, in the fitness space, 
of the distance between the infeasible solution with the lower objective and the 
individual in the population with the highest objective function.

The aim of the technique is not to reject completely individuals that 
are not totally feasible, or in other words that slightly violate some
of the constraints. The actual implementation of the self-adaptive 
method in PaGMO/PyGMO is based on population fitness evaluation 
modification that adapts at each algorithm generation.
The main advantages of this techniques are that it does not require parameter 
tuning and can be used also without a feasible solution in the initial
population. 

Application
###########
The problem considered here is the problem g01 from the Congress on 
Evolutionary Computation 2006 (CEC2006). This problem has a quadratic
objective function with nine linear inequality constraints. 
In the optimum six constraints are active. 

To solve it, a simple genetic algorithm with Gray encoding has been used
in the meta-algorithm defined by the constraints handling technique described above. 
The meta-algorithm wraps the original 
algorithm that was only dedicated to solve unconstrained problems. It assignes to the population
a modified problem that performes the adaptive penalization. In the example 70 individuals 
are considered into the population. 

Step by step, the PyGMO library is imported and the
populations size and the number of generation for the algorithm are initialized.

.. code-block:: python

   In [1]: from PyGMO import *
   In [2]: pop_size = 70
   In [3]: n_gen = 20000

Then the problem is created, a corresponding population randomly initialized and the genetic algorithm instantiated.

.. code-block:: python

   In [4]: prob = problem.cec2006(1)
   In [5]: pop = population(prob,pop_size)
   In [6]: algo = algorithm.sga_gray(gen=1,cr=0.9,m=0.004,elitism=0,mutation=algorithm.sga_gray.mutation.UNIFORM,selection=algorithm.sga_gray.selection.ROULETTE,crossover=algorithm.sga_gray.crossover.SINGLE_POINT)

In this example the genetic algorithm has a crossover probability of 90%, a mutation probability 
of 0.4%. Furthermore, it uses a uniform mutation, a roulette selection 
and a single point crossover. For more details on this algorithm, please
refer to its documentation.

The important point here is the number of generation of the algorithm.
In fact, as the number of generation is governed by the meta-algorithm
and that the meta-algorithm modifies the fitness of the population at
each iterations, the number of generation of the algorithm must be set
to 1 to work properly. The same comment holds for the elitism.

The next step consists in creating the meta-algorithm that implements the self-adaptive constraints handling technique presented.

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

   In [12]: print(prob.best_x)
   In [13]: print(prob.best_f)
   In [14]: print(prob.best_c)
   Out [2]:
   ((1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0, 3.0, 3.0, 1.0),)
   ((-15.0,),)
   ((0.0, 0.0, 0.0, -5.0, -5.0, -5.0, 0.0, 0.0, 0.0),)

Even if the solutions are really close to the optimum, the exact 
same performances of this algorithm as described by R. Farmani 
and J. A. Wright in their paper could not be retrieved with our 
configuration. This might be due to a different implementation of the heuristic technique used for the optimization.
