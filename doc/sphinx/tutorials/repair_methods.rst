.. _repair_methods:

=======================================================================
Repair methods
=======================================================================
Repair methods is a constraints handling technique that consists in 
repairing infeasible individuals of the population to make 
them approaching the feasible region. In this tutorial, we will learn how to solve a 
constrained problem with this technique.

Method
##########
The repairing method implemented in PaGMO/PyGMO is extremely simple.
It is based on a meta-algorithm where to the initial problem is 
associated an unconstrained problem obtained simply removing the 
constraints an leaving unvaried the original objective function. 
The problem is solved for optimality and after a predefined number of generations, 
infeasible individuals are repaired. The name of this meta-algorithm
is CORE (Constrained Optimization by Random Evolution). 
The repairing method uses a simple gradient descent methods
to repri the infeasibility of the individuals. Thus, the
meta algorithm takes two algorithms, one to evolve the main population towards optimality,
and one for repairing. 

Application
###########
The problem considered here is the problem g05 from the Congress on 
Evolutionary Computation 2006 (CEC2006). This problem has a cubic
ojective function with two linear inequality and three non linear equality 
constraints. To solve this problem the Differential
Evolution (DE) is chosen as the main evolutionary technique and a simplex method as the repair 
algorithm. The simplex method used is included into the GSL library. This means that PyGMO/PaGMO
needs to be compiled with the GSL option activated. The population 
contains 90 individuals.

First import the PyGMO library and choose the populations size and the
number of generation for the meta-algorithm.

.. code-block:: python
   
   In [1]: from PyGMO import *
   In [2]: pop_size = 90
   In [3]: n_gen = 1000
   In [4]: n_repair_gen = 100

Then instantiate the algorithms you wish to use. The generation of the
first main algorithm must be set to 1 as the number of iterations is
driven by the meta-algorithm.

.. code-block:: python

   In [5]: algo_1 = algorithm.de(gen = 1, xtol=1e-30, ftol=1e-30)
   In [6]: algo_repair = algorithm.gsl_nm2(max_iter = n_repair_gen, step_size = 0.02, tol = 1e-8)

Select the problem and associate a population to this problem.

.. code-block:: python

   In [7]: prob = problem.cec2006(5)
   In [8]: pop = population(prob,pop_size)

Creates the meta-algorithm with these informations.

.. code-block:: python

   In [9]: algo_meta = algorithm.cstrs_core(algorithm = algo_1, repair_algorithm = algo_repair, gen = n_gen, repair_frequency = 10, repair_ratio = 1., f_tol = 1e-15, x_tol = 1e-15)

A repairing frequency of 10 generations is selected, and 
all the infeasible individuals are repaired, hence the ratio is set to 1.

The evolution is then performed.

.. code-block:: python

   In [10]: pop = algo_meta.evolve(pop)

And finally, print the solutions.

.. code-block:: python

   In [11]: print(pop.champion.x)
   In [12]: print(pop.champion.f)
   In [13]: print(pop.champion.c)
   Out [1]:
   (746.9608574722797, 955.4049402203846, 0.07150528169781202, -0.4189587300214994)
   (5149.854514758014,)
   (6.983932576076768e-05, -3.644692628768098e-05, 2.23429769903305e-05, -0.05953598828068862, -1.0404640117193114)

As a comparison, you can print the best known solution for this
particular problem:

.. code-block:: python

   In [11]: print(prob.best_x)
   In [12]: print(prob.best_f)
   In [13]: print(prob.best_c)
   Out [2]:
   ((679.9451482970287, 1026.066976000047, 0.11887636909441043, -0.39623348521517826),)
   ((5126.4967140071,),)
   ((9.999999997489795e-05, 9.999999997489795e-05, 9.999999997489795e-05, -0.03489014569041138, -1.0651098543095887),)

Note that you might need to multiple run this tutorial to get a
feasible solution.

If for any reason you wish to repair by hand any individual of the 
population, you can proceed as follow:

.. code-block:: python

   In [11]: pop.repair(0,algo_repair)

In this case, the individual at index 0 is repaired with the algorithm 
algo_repair.
