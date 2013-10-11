.. _multi_objective_transformation:

=======================================================================
Multi-objective transformation
=======================================================================
The multi-objective transformation is a technique to solve constrained
problem by transforming the initial single (or multi) objective constrained problem 
into a multi-objective unconstrained problem with the number of objectives
equal to the sum of the number of objectives and the number of constraints. 
In this tutorial, we are going to solve a single objective constrained problem
using this technique.

Method
##########
The multi-objective transformation consists in treating each
constraints of a constrained optimization problem as objectives of a multi-objective
problem, the first objective(s) being the cost function(s) itself. Three
different implementation of this method are available in PaGMO/PyGMO.

The first one called OBJ_CSTRS splits all the constraints into 
different objectives, generating a problem with *nobj+m* objectives (where
*nobj* is the number of objectives and *m* is the number of constraints). 
The transformation adopted is the one of Coello Coello in the CHVEGA algorithm where the *nobj+j* objective 
is the j-th constraints violation (if the j-th constraint is violated), 
otherwise it is the overall number of violated constraints, otherwise if the solution 
is feasible it is the aggregation of the objectives. This formulation makes sense if the 
multiobjective algorithm is based on a population decomposition according to the different 
objectives as VEGA or PADE.
The second one called OBJ_CSTRSVIO has *nobj+1* objectives, the first nobj 
are the objective functions to be optimized and the
last objective is the sum of the constraints violations. 
The third methodology has *nobj+2*
objectives, the first nobj are the objective functions, the second is the
sum of inequality constraints violation and the third one is the sum of equality
constraints violation. 

In all cases, when a constraint is satisfied, the 
associated objective is either a measure of the number of the violated
constraints, if at least one constraint is violated, or the objective 
function by itself. That way, each objectives commonly helps to reduce
the number of constraints violation to find feasible solutions.

Application
###########
The problem considered here is the problem g04 from the Congress on 
Evolutionary Computation 2006 (CEC2006). This problem has a quadratic
objective function with six non linear inequality constraints. 
In the optimum two constraints are active.

To solve it, in PaGMO/PyGMO, the
meta-problem that modifies the original constrained problem into a 
multi-objective one is used. Then a multi-objective optimizer called 
VEGA is used to solve it. This optimizer is a simple extention of simple genetic algorithm 
to solve multi-objective algorithms. The population considered is composed 
by 280 individulas (40 individuals per objectives: (6+1) * 40 = 280 since VEGA equally
divides the original population into a number of subpopulations equal to the number of objectives). 

Step by step, the PyGMO library is imported and the
populations size and the number of generation for the algorithm are initialized.

.. code-block:: python

   In [1]: from PyGMO import *
   In [2]: pop_size = 280
   In [3]: n_gen = 500

Then the problem, the meta-problem (using the first methodology OBJ_CSTRS) and the algorithm are created.

.. code-block:: python

   In [4]: prob = problem.cec2006(4)
   In [5]: prob_mo = problem.con2mo(prob,problem.con2mo.method.OBJ_CSTRS)
   In [6]: algo_mo = algorithm.vega(gen = n_gen)

The next step consists in creating the population and evolving it.

.. code-block:: python

   In [7]: pop = population(prob_mo, pop_size)
   In [8]: pop = algo_mo.evolve(pop)

From that point, the obtained solutions can be printed.

.. code-block:: python

   In [9]: print(pop.champion.x)
   In [10]: print(prob.objfun(pop.champion.x))
   In [11]: print(prob.compute_constraints(pop.champion.x))
   Out [1]:
   (78.09947909827132, 33.30065259242254, 30.894527032856043, 44.013876005779935, 35.07762608948753)
   (-30476.229905572058,)
   (-0.2613258099525382, -91.73867419004746, -11.284338550618045, -8.715661449381955, -4.980246388885227, -0.019753611114772696)

The solution found by this algorithm is close to the optimum value:

.. code-block:: python

   In [12: print(prob.best_x)
   In [13]: print(prob.best_f)
   In [14]: print(prob.best_c)
   Out [2]:
   ((78.0, 33.0, 29.9952560256816, 45.0, 36.77581290578821),)
   ((-30665.538671783317,),)
   ((0.0, -92.0, -11.159499691073137, -8.840500308926863, -4.9999999999999964, -3.552713678800501e-15),)

Due to the stochastic behavior of the algorithm in a single run it might not converge to the exact optimum. We invite you
to play with the population size and the number of generations 
to see the behavior of this constraints handling technique.
