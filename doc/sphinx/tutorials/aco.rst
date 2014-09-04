.. _aco:

=======================================================================
Solving Traveling Salesman Problems
=======================================================================

Given a list of cities and their pairwise distances, the task is to find the 
shortest possible tour that visits each city exactly once.

TSP can be modeled as a graph, such that cities are the graph's vertices, 
paths are the graph's edges, and a path's distance is the edge's weight. 
A TSP tour becomes a Hamiltonian cycle, and the optimal TSP tour is the 
shortest Hamiltonian cycle.

In the symmetric TSP, the distance between two cities is the same in 
each opposite direction, forming an undirected graph. 
This symmetry halves the number of possible solutions. 
In the asymmetric TSP, paths may not exist in both directions 
or the distances might be different, forming a directed graph.
Traveling salesman problems are formulated in PyGMO as integer linear 
programming problems.
Each vertex is labeled with an index number, where n is the total number of vertices.
:math:`X_{i,j}` is the square binary matrix, with 

.. math::

        X_{i,j} = 
        \begin{cases}
            1, & \text{if path goes from city i to city j}\\
            0, & \text{otherwise}
        \end{cases}

:math:`C_{i,j}` is the square adjacency matrix containing the distances (weights)
between city i and city j, with i and j :math:`\in = \{1, 2, .., n\}`. 
If the problem is symmetric, then the matrix is symmetric, meaning it's equal to it's transpose.
For :math:`i = \{1, 2, .., n\} \quad u_i` is an artificial variable which has to satisfy the inequality constraints.
Then, the TSP can be written as follows:

.. math::

        \min \sum_{i=0, i \neq j}^n \sum_{j=0}^n C_{i,j} X_{i,j}

        0 \le X_{i,j} \le 1 \quad i, j = \{1, 2, .., n\}

        u_i \in \mathbb{Z} \quad i = \{1, 2, .., n\}

        \sum_{i=0, i \neq j}^n X_{i,j} = 1 \quad j = \{1, 2, .., n\}
        
        \sum_{j=0, i \neq j}^n X_{i,j} = 1 \quad i = \{1, 2, .., n\}

        u_i - u_j + n * X_{i,j} \le n - 1 \quad i \le i \neq j \le n

        
For more information visit: http://en.wikipedia.org/wiki/Travelling_salesman_problem#Integer_linear_programming_formulation

Instantiating a TSP problem
###########################

Here we give two examples on how to instantiate a TSP problem from a weights matrix
or a TSPLIB XML file http://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/

From weights
------------

The weights matrix is a two dimensional vector of doubles, e.g.

.. code-block:: python

        weights = [[0, 1, 2], [3, 0, 5], [6, 7, 0]]

PyGMO ignores main diagonal elements. Notice in the example that the diagonal elements are zero.

In the first example we create a list of list (matrix) 
and print the instantiated tsp problem to console:

.. code-block:: python

        from PyGMO.problem import tsp
        
        # creating a list of lists
        list2D = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
        
        # creating an instance of a TSP problem
        tsp_instance = tsp(list2D)
        
        # printing the TSP instance to console
        print tsp_instance

And the resulting TSP problem looks like this. Again, notice there are no
connections between a vertex and itself, regardless of what we put in the main diagonal,
in our case the values 0, 4 and 8 are ignored.

.. code-block:: python

        Problem name: Traveling Salesman Problem
	        Global dimension:			6
	        Integer dimension:			6
	        Fitness dimension:			1
	        Constraints dimension:			8
	        Inequality constraints dimension:	2
	        Lower bounds: [0, 0, 0]
	        Upper bounds: [1, 1, 1]
	        Constraints tolerance: []
        The Boost Graph (Adjacency List): 
        Vertices = { 0 1 2 }
        Edges (Source, Target) = Weight : 
        (0, 1) = 1
        (0, 2) = 2
        (1, 0) = 3
        (1, 2) = 5
        (2, 0) = 6
        (2, 1) = 7
        
Now, if you didn't read the wikipedia page and were wondering about the numbers, here's whey they come from:

Global dimension :math:`= n*n-n = n*(n-1) = 3 * 2 = 6`

Integer dimension :math:`= n*(n-1) = 3 * 2 = 6`

Fitness dimension (max) :math:`= 1 \in [0,1]`

Constraints dimension (global) = Equality + Inequality :math:`= 2*n + (n-1)*(n-2) = n*(n-1)+2 = 3 * 2 + 2 = 8`

Equality constraints dimension :math:`= 2*n`

Inequality constraints dimension :math:`= (n-1)*(n-2) = 2 * 1 = 2`


From TSPLIB XML
---------------

In this second example we will be loading an TSPLIB XML file from the current folder (pwd in linux).

In order to load an XML file, we use the utility function PyGMO.util.tsp.read_tsplib('file.xml')
which returns a weights matrix (list of list) such as the one defined above.

.. code-block:: python

        from PyGMO.util import tsp as tsputil
        from PyGMO.problem import tsp

        # importing the XML file
        weights = tsputil.read_tsplib('burma14.xml')

        # printing the weights matrix
        tsputil.print_matrix(weights)

        # creating a tsp problem from the imported weights matrix
        tsp_instance = tsp(weights)

        # printing the tsp problem details to console
        print tsp_instance

The imported matrix, notice it is symmetric.

.. code-block:: python

        [[    0.   153.   510.   706.   966.   581.   455.    70.   160.   372.   157.   567.   342.   398.]
         [  153.     0.   422.   664.   997.   598.   507.   197.   311.   479.   310.   581.   417.   376.]
         [  510.   422.     0.   289.   744.   390.   437.   491.   645.   880.   618.   374.   455.   211.]
         [  706.   664.   289.     0.   491.   265.   410.   664.   804.  1070.   768.   259.   499.   310.]
         [  966.   997.   744.   491.     0.   400.   514.   902.   990.  1261.   947.   418.   635.   636.]
         [  581.   598.   390.   265.   400.     0.   168.   522.   634.   910.   593.    19.   284.   239.]
         [  455.   507.   437.   410.   514.   168.     0.   389.   482.   757.   439.   163.   124.   232.]
         [   70.   197.   491.   664.   902.   522.   389.     0.   154.   406.   133.   508.   273.   355.]
         [  160.   311.   645.   804.   990.   634.   482.   154.     0.   276.    43.   623.   358.   498.]
         [  372.   479.   880.  1070.  1261.   910.   757.   406.   276.     0.   318.   898.   633.   761.]
         [  157.   310.   618.   768.   947.   593.   439.   133.    43.   318.     0.   582.   315.   464.]
         [  567.   581.   374.   259.   418.    19.   163.   508.   623.   898.   582.     0.   275.   221.]
         [  342.   417.   455.   499.   635.   284.   124.   273.   358.   633.   315.   275.     0.   247.]
         [  398.   376.   211.   310.   636.   239.   232.   355.   498.   761.   464.   221.   247.     0.]]
        
And finally, the output for printing the TSP problem instance:

.. code-block:: python

    Problem name: Traveling Salesman Problem
	Global dimension:			182
	Integer dimension:			182
	Fitness dimension:			1
	Constraints dimension:			184
	Inequality constraints dimension:	156
	Lower bounds: [0, 0, ..., 0]
	Upper bounds: [1, 1, ..., 1]
	Constraints tolerance: [0, 0, ..., 0]

        The Boost Graph (Adjacency List): 
        Vertices = { 0 1 2 3 4 5 6 7 8 9 10 11 12 13 }
        Edges (Source, Target) = Weight : 
        (0, 1) = 153.0
        # [..snip..]
        (13, 12) = 247.0
        
Solving Using Ant Colony Optimization
#####################################

Here we solve the 'kroA100.xml' problem file from TSPLIB,
using three variants of Ant Colony Optimization algorithms.

Simple Ant Colony
-----------------

In Ant Colony Optimization algorithms, ants are randomly placed in randomly 
starting positions at each cycle start. As each ant travels, it deposits pheromone.
Ants decide which path to take according to the distance between vertices (edges)
and the amount of pheromone deposited on each edge. Initially this is set to a small constant.

The balance between distance and pheromone is controlled by two parameters, alpha and beta.
When an ant makes a decision between the current edge and the next, the probability of the transition is:
p(next|current) = ( pheromone(next, current) ^ alpha * 1/distance(next, current) ^ beta) ) / Sum of all possibilities

Since shorter paths are traversed more often, more pheromone is deposited on them,
thus on the long term this greedy strategy converges to an optimum.
To allow the algorithm to explore the search space as much as possible, 
we force the pheromone to evaporate, as time passes, using a constant rho.

We run the algorithm for 1000 cycles to allow it to converge.
The number of ants is usually set equal to the number of cities, in our case 100.

.. code-block:: python

    from PyGMO import *
    from PyGMO.util import tsp as tsputil

    # Importing tsp file
    xml = tsputil.read_tsplib('kroA100.xml')

    # Instantiatng a tsp problem
    prob = problem.tsp(xml)

    # Creating population of 99 individuals.
    # The number of ants must be greater than the population.
    pop = population(prob, 99)

    # ---- Instantiating a Simple Ant Colony System ----
    # We run the algorithm for 1000 cycles to allow it to converge.
    # This parameter is the only one required, the others are optional.
    # The complete constructor signature is: aco(no_cycles, no_ants, rho, alpha, beta)
    #
    # The number of ants is usually set equal to the number of cities, in our case 100.
    # Rho, the third argument, controls how fast the pheromone evaporates as time passes.

    algo_simple = algorithm.aco(1000, 100, 0.5)

    # Aditionally, all these properties have getters and setters.
    print algo_simple.gen
    print algo_simple.ants
    print algo_simple.rho

    # We haven't talked yet about alpha and beta.
    # These two parameters control the way that ants make local decisions,
    # as a compromise between the length of an edge versus the amount of pheromone deposited on it.
    # When an ant makes a decision between taking edge A or B, the probability of the transition is:
    # p(next) = ( pheromone(A, B) ^ alpha * 1/distance(A,B) ^ beta) ) / Sum of all possibilities
    #
    # Let's explicitly call the setters for these two.
    # Good values are usually 1 for alpha and between 2 and 5 for beta
    algo_simple.alpha = 1
    algo_simple.beta = 2

    # To actually run the algorithm, we call the evolve method
    result_simple = algo_simple.evolve(pop)

    # The champion contains the fitness, which is the length of the shortest tour
    print result_simple.champion.f


Elite Ant Colony
-----------------

In this variant, the best ant deposits additional pheromone in each cycle.
This amount is set as parameter 'e' which is usually set equal to the number of ants.
Since the best ant (shortest distance traveled) deposits much more pheromone, 
it increases the likelihood that other ants will take the same path on the next cycles. 
Remember alpha and beta from Simple ACO. These along with the initialization of the
pheromone influence how fast and accurate the algorithm will converge.


.. code-block:: python

    from PyGMO import *
    from PyGMO.util import tsp as tsputil

    # Importing tsp file
    xml = tsputil.read_tsplib('kroA100.xml')

    # Instantiatng a tsp problem
    prob = problem.tsp(xml)

    # Creating population of 99 individuals.
    # The number of ants must be greater than the population.
    pop = population(prob, 99)

    # ---- Instantiating an Elite Ant Colony System
    # In this variant, the best ant for each cycle deposits additional pheromone.
    # This amount is set as parameter 'e' which is usually set equal to the number of ants.
    algo_elite = algorithm.aco_elite(1000, 100, 0.5, 100)

    # The elite constant also has setters and getters
    print algo_elite.e

    # We run the evolve method
    result_elite = algo_elite.evolve(pop)

    # The champion contains the fitness, which is the length of the shortest tour
    # Since the best ant deposits much more pheromone, it increases the likelihood that
    # other ants will take the same path on next runs. Remember alpha and beta from Simple ACO.
    # We print the whole details for the champion. 
    # Note that the tour is in binary (chromosome) format.
    print result_elite.champion

Rank-Based Ant Colony
---------------------

In rank-based ACO the idea is similar to the Elite ACO. In addition,
each ant deposits an amount of pheromone which is proportional to it's rank.
The shorter an ant makes a tour, the better it's rank, the more pheromone it deposits.

Typically this form of algorithm has a constant integer which controls the number of ranks. 
In our implementation, the number of ranks is equal to the number of unique solutions (tour lengths found) in each cycle.
In this case, the elite constant 'e' controls how much pheromone each ant deposits.
The best ant(s) deposit 100% of 'e' while the worst does not deposit any.
In other words, the best ants behave just like in the elite ant colony,
while the other ants deposit according to this formula: 
percentage(ant) = (no_ranks - tour_distance(ant) )/no_ranks
Where no_ranks is the number of unique solutions found in a cycle.


.. code-block:: python

    from PyGMO import *
    from PyGMO.util import tsp as tsputil

    # Importing tsp file
    xml = tsputil.read_tsplib('kroA100.xml')

    # Instantiatng a tsp problem
    prob = problem.tsp(xml)

    # Creating population of 99 individuals.
    # The number of ants must be greater than the population.
    pop = population(prob, 99)

    # ---- Instantiating an Rank-based Ant Colony System
    # Each ant deposits an amount of pheromone which is proportional
    # to it's rank (length of tour) that it has performed.
    # In other words, the shorter the tour, the more pheromone an ant deposits.
    # The elite constant 'e' controls how much pheromone each ant deposits.
    # The best ant(s) deposit 100% 'e' while the worst does not deposit any.
    # We need to specify the number of cycles in order to instantiate any ACO algorithm.

    algo_rank = algorithm.aco_rank(1000)

    # Let's set the parameters using the setters in this case.
    algo_rank.ants = 100
    algo_rank.rho = 0.3
    algo_rank.e = 150

    # Running the evolve method
    result_rank = algo_rank.evolve(pop)

    # Print the results
    print result_rank.champion.f

