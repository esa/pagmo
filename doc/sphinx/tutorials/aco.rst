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
connections between a vertice and itself, regardless of what we put in the main diagonal,
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


