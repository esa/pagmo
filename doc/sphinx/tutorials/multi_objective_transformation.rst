.. _multi_objective_transformation:

=======================================================================
Multi-objective transformation
=======================================================================
The multi-objective transformation is a technique to solve constrained
problem by transforming the initial problem into a multi-objective
problem. In this tutorial, we are going to solve a constrained problem
using this technique.

Method
##########
The multi-objective transformation consists in treating each
constraints of a constrained problem as objectives of a multi-objective
problem, the first objective being the cost function by itself. Three
different implementation of this method are available in PaGMO/PyGMO.
The first one called OBJ_CSTRS splits all the constraints into 
different objectives. The second one called OBJ_CSTRSVIO has two
objectives, the first one is the function to optimize and the
second objective is the sum of the constraints. The third one has three
objectives, the first one is the objective function, the second the
sum of inequality constraints and the third one the sum of equality
constraints. In all cases, when a constraint is satisfied, the 
associated objective is either a measure of the number of the violated
constraints, if at least one constraint is violated, or the objective 
function by itself. That way, each objectives commonly helps to reduce
the number of constraints violation to find feasible solutions.

Application
###########
We consider here the constrained problem g04 from the Congress on 
Evolutionary Computation 2006 (CEC2006). This problem is a quadratic
problem with 6 nonlinear inequality constraints and 2 active
constraints. To solve it, we use a multi-objective optimizer called 
VEGA. This optimizer is a simple extention of simple genetic algorithm 
to solve multi-objective algorithms. It is coupled with this 
constraints handling technique. To do so in PaGMO/PyGMO, we use a 
meta-problem that modifies the original constrained problem into a 
multi-objective one. For the population, we consider 40 individuals 
per objectives or a total of 280 individuals: (6+1) * 40 = 280. In 
fact, this constraints handling technique automatically splits the 
initial population to equaly split it in the number of objectives.

Step by step, we first import the PyGMO library and choose the
populations size and the number of generation for the algorithm.
.. code-block:: python
   In [1]: from PyGMO import *
   In [2]: pop_size = 280
   In [3]: n_gen = 500

Then we create the problem, the meta-problem and the algorithm.
.. code-block:: python
   In [4]: prob = problem.cec2006(4)
   In [5]: prob_mo = problem.con2mo(prob,problem.con2mo.method.OBJ_CSTRS)
   In [6]: algo_mo = algorithm.vega(gen = n_gen)

The next step consists in creating the population and evolving it.
.. code-block:: python
   In [7]: pop = population(prob_mo, pop_size)
   In [8]: pop = algo_mo.evolve(pop)

From that point, we can print the solutions.
.. code-block:: python
   In [9]: print(pop.champion.x)
   In [10]: print(pop.champion.f)
   In [11]: print(pop.champion.c)

   Out [1]:
(78.09947909827132, 33.30065259242254, 30.894527032856043, 44.013876005779935, 35.07762608948753)
(-30476.229905572058,)
(-0.2613258099525382, -91.73867419004746, -11.284338550618045, -8.715661449381955, -4.980246388885227, -0.019753611114772696)

The solution found by this algorithm is close to the optimum given
with the following:
.. code-block:: python
   In [12: print(prob.best_x)
   In [13]: print(prob.best_f)
   In [14]: print(prob.best_c)

   Out [2]:
((78.0, 33.0, 29.9952560256816, 45.0, 36.77581290578821),)
((-30665.538671783317,),)
((0.0, -92.0, -11.159499691073137, -8.840500308926863, -4.9999999999999964, -3.552713678800501e-15),)

Of course, you will not get exactly the same results as we did 
due to the stochastic behavior of the algorithm. We invite you
to play with the population size and the number of generations 
to see the behavior of this constraints handling technique.
