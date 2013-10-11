.. _death_penalty:

=======================================================================
Death penalty
=======================================================================

Death penalty is a static penalty constraints handling technique:
individuals that are not satisfying the constraints are penalized with 
a huge constraint value. The technique assumes minimization and it can be 
applyed to both single and multi-objective constrained problems.
In this tutorial we will learn how to use this specific constraints 
handling technique using PaGMO/PyGMO.

Method
##########
The death penalty technique is implemented through a meta-problem. 
In our case, the meta-problem takes a constrained problem, removes 
its constraints and penalizes the original objective function value 
with a factor that is a measure of the point infeasibility.
The derived meta-problem can then be solved by any optimization 
algorithm for unconstrained optimization. 

In PaGMO/PyGMO, two different death penalty techniques are implemented.
The first one is the simple death penalty method where infeasibility is
penalized assigning to the objective value the maximum machine-finite value 
[boost::numeric::bounds<double>::highest()].
The second one is the Kuri technique where the same penalization value 
is applyed according to a rate of constraints satisfaction.

Application
###########
The problem considered here is the problem g04 from the Congress on 
Evolutionary Computation 2006 (CEC2006). This problem has a quadratic
objective function with six non linear inequality constraints. 
In the optimum two constraints are active. 

To solve this problem two heuristic techniques are used: the Differential
Evolution (DE) and the Simple Genetic Algorithm (SGA). The population size
is 70 and 5000 generations are performed, which gives a total number of 
functions evaluation equal to 350000. These parameters can be set as preferred in
the constructor of problem and algorithms. 
To get statistical meaning 25 runs are performed per each algorithm and mean objective value,
the best optimum found and the standard deviation are computed.

The code is very explicit by itself. Copy it to a file named run.py
in example, and run this file with python.

.. code-block:: python

   from PyGMO import *
   from numpy import *

   n_trials = 25
   pop_size = 70
   n_gen = 5000

   prob_cec = problem.cec2006(4)

   prob_list = []
   prob_list.append(problem.death_penalty(prob_cec,problem.death_penalty.method.SIMPLE))
   prob_list.append(problem.death_penalty(prob_cec,problem.death_penalty.method.KURI))

   print('\n-----------------------------')
   print('\nOriginal problem: ')
   print(prob_cec)
   print('\nMeta-problem 0: ')
   print(prob_list[0])
   print('\n-----------------------------')
   print('\nBest decision vector: ' +  str(prob_cec.best_x[0]) )
   print('\nBest fitness: ' +  str(prob_cec.best_f[0]) )
   print('\n-----------------------------')

   algo_list = [algorithm.de(gen = n_gen, xtol=1e-30, ftol=1e-30), algorithm.sga(gen = n_gen)]

   print('\nTrials: ' + str(n_trials) + ' - Population size: ' + str(pop_size) + ' - Generations: ' + str(n_gen))

   for prob in prob_list:
       print('\n-----------------------------')
       print('\nTesting problem: ' + prob.get_name() + ' with population size: ' +  str(pop_size) )
            
       for algo in algo_list:
           print('\n' + str(algo))
           best = []
           best_x = []
           for i in range(0,n_trials):
               isl = island(algo,prob,pop_size)
               isl.evolve(1)
               isl.join()
               if prob_cec.feasibility_x(isl.population.champion.x):
                   best.append(isl.population.champion.f)
                   best_x.append(isl.population.champion.x)

           print(' Best:\t' + str(min(best)[0]))
           print(' Mean:\t' + str(mean(best)))
           print(' Std:\t' + str(std(best)))

           print(' Example of found x:\t' + str(best_x[0]))
           print(' Example of found constraints:\t' + str(prob_cec.compute_constraints((best_x[0]))))

In this example, we have used two different meta-problems one for each death
penalty techniques, both contained in the prob_list variable. By looking at
the output given by the original problem and the first meta problem, 
it is important to see that the original constrained problem output has 6 
constraints while the meta-problem doens't have any as expected.

If run directly into python, you would get the following output:

.. code-block:: python

   Out[1]:
   -----------------------------

   Original problem: 
   Problem name: CEC2006 - g4
   	Global dimension:			5
   	Integer dimension:			0
   	Fitness dimension:			1
   	Constraints dimension:			6
   	Inequality constraints dimension:	6
   	Lower bounds: [78, 33, 27, 27, 27]
   	Upper bounds: [102, 45, 45, 45, 45]
   	Constraints tolerance:			[0, 0, 0, 0, 0, 0]


   Meta-problem: 
   Problem name: CEC2006 - g4 [death_penalty, method_SIMPLE ]
   	Global dimension:			5
   	Integer dimension:			0
   	Fitness dimension:			1
   	Constraints dimension:			0
   	Inequality constraints dimension:	0
   	Lower bounds: [78, 33, 27, 27, 27]
   	Upper bounds: [102, 45, 45, 45, 45]
   	Constraints tolerance:			[]


   	Constraints handled with death penalty, method SIMPLE 

   -----------------------------
   
   Best decision vector: (78.0, 33.0, 29.9952560256816, 45.0, 36.77581290578821)

   Best fitness: (-30665.538671783317,)

   -----------------------------
   Trials: 25 - Population size: 70 - Generations: 3000
   -----------------------------

   Testing problem: CEC2006 - g4 [death_penalty, method_SIMPLE ] with population size: 70

   Algorithm name: Differential Evolution - gen:3000 F: 0.8 CR: 0.9 variant:2 ftol:1e-30 xtol:1e-30
    Best:	-30665.5386718
    Mean:	-30665.5386718
    Std:	3.8500748794e-12
    Example of found x:	(78.0, 33.0, 29.99525602568156, 44.99999999999992, 36.775812905788335)
    Example of found constraints:	(0.0, -92.0, -11.159499691073108, -8.840500308926892, -5.0, 0.0)

   Algorithm name: A Simple Genetic Algorithm - gen:3000 CR:0.95 M:0.02 elitism:1 mutation:GAUSSIAN (0.1) selection:ROULETTE crossover:EXPONENTIAL 
    Best:	-30645.1250077
    Mean:	-30602.5590724
    Std:	74.4327472986
    Example of found x:	(78.01739822664655, 33.14346829880361, 30.306842224157915, 44.81436783797341, 36.06696016578172)
    Example of found constraints:	(-0.09005427354273365, -91.90994572645727, -11.213189038009048, -8.786810961990952, -4.999961858160862, -3.814183913775082e-05)

   -----------------------------

   Testing problem: CEC2006 - g4 [death_penalty, method_KURI ] with population size: 70

   Algorithm name: Differential Evolution - gen:3000 F: 0.8 CR: 0.9 variant:2 ftol:1e-30 xtol:1e-30
    Best:	-30665.5386718
    Mean:	-30665.5386718
    Std:	4.05108175097e-12
    Example of found x:	(78.0, 33.0, 29.99525602568155, 44.99999999999988, 36.77581290578837)
    Example of found constraints:	(0.0, -92.0, -11.15949969107308, -8.84050030892692, -5.0, 0.0)

   Algorithm name: A Simple Genetic Algorithm - gen:3000 CR:0.95 M:0.02 elitism:1 mutation:GAUSSIAN (0.1) selection:ROULETTE crossover:EXPONENTIAL 
    Best:	-30659.9680031
    Mean:	-30616.4063446
    Std:	28.812620113
    Example of found x:	(78.04964811877262, 33.07572602463478, 30.150381195899435, 44.965160220602215, 36.46285810358295)
    Example of found constraints:	(-0.03509668936189314, -91.9649033106381, -11.170496308902116, -8.829503691097884, -4.989167173895257, -0.010832826104742566)


This example shows that the 
Differential Algorithm with both the simple and the Kuri death
penalty methods converges to the best known decision vector.
Furthermore, it is noticeable that on this case, the Differential
Evolution algorithm outperforms the Simple Genetic Algorithm.

Please note that, from the literature, it is stated that the death penalty constraints handling method is
not well suited to solve highly constrained problems or problems
with equality constraints. Indeed if at initialization no feasible individuals 
are present, all the individuals will have assigned the same penalty value and the 
algorithm will have no preferible direction to follow and will proceeds just with random tentatives.
