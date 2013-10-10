.. _decomposition:

================================================================
Decomposition
================================================================

In this tutorial we will learn how to use the decompose meta-problems to solve a multi-objective problem.
The decompose meta-problem transforms a multi-objective problem into a single-objective one having as
fitness function a convex combination (defined by a weight vector) of the original objectives. 
Let us start creating a decomposed problem from a multi-objective one.

.. code-block:: python
   
	In [1]: from PyGMO import *
	In [2]: orig_prob = problem.zdt(1,10)
	In [3]: prob = problem.decompose(problem = orig_prob, weights = [0.5, 0.5])

In this way the 2 objectives of the original problem are equally weighted. If we don't define the
weight vector, then it is randomly generated.

We then proceed by solving the new decomposed problem using a single-objective optimization algorithm.

.. code-block:: python

	In [4] alg = algorithm.jde(50)
	In [5] pop = population(prob, 200)
	In [6] for i in xrange(5):
		    pop = alg.evolve(pop)
		    print "Generation ", i
		    print "Distance from Pareto Front (p-distance): " , orig_prob.p_distance(pop) 
		    print "Original fitness: " , orig_prob.objfun(pop.champion.x)
		    print "Decomposed fitness: " , pop.champion.f

We see how the fitness of the new problem is equal to the average of the original 2 objectives 

.. code-block:: python

	Out [7] Distance from Pareto Front (p-distance):  0.569852362819
		Original fitness:  (0.376410003683959, 0.6055773596287406)
		Decomposed fitness:  (0.4909936816563498,)

		Generation  1
		Distance from Pareto Front (p-distance):  0.0962829285609
		Original fitness:  (0.25292383041734323, 0.528506927611418)
		Decomposed fitness:  (0.39071537901438064,)

		Generation  2
		Distance from Pareto Front (p-distance):  0.0164014486696
		Original fitness:  (0.25426237266331964, 0.499545293860021)
		Decomposed fitness:  (0.3769038332616703,)

		Generation  3
		Distance from Pareto Front (p-distance):  0.00251788007034
		Original fitness:  (0.24129082628054865, 0.5095238228445969)
		Decomposed fitness:  (0.37540732456257275,)

		Generation  4
		Distance from Pareto Front (p-distance):  0.000420874501332
		Original fitness:  (0.2500151198152518, 0.5000990342618756)
		Decomposed fitness:  (0.3750570770385637,)

		Generation  5
		Distance from Pareto Front (p-distance):  7.23681824182e-05
		Original fitness:  (0.24995048982014567, 0.5000683842135499)
		Decomposed fitness:  (0.3750094370168478,)

