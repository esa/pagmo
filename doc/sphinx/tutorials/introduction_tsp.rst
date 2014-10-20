.. _introduction_tsp:

=======================================================================
Introduction
=======================================================================

Given a list of cities and their pairwise distances, the task is to find the 
shortest possible tour that visits each city exactly once.

In the symmetric TSP, the distance between two cities is the same in 
each opposite direction, forming an undirected graph. 
This symmetry halves the number of possible solutions. 
In the asymmetric TSP, paths may not exist in both directions 
or the distances might be different, forming a directed graph.

Traveling salesman problems can be formulated (encoded) in PyGMO in three different ways.

TSP as  integer linear programs (FULL encoding)
-----------------------------------------------

In this encoding each edge is associated to a variable :math:`X_{i,j}`

:math:`X_{i,j}` is the square binary matrix, with 

.. math::

        X_{i,j} = 
        \begin{cases}
            1, & \text{if path goes from city i to city j}\\
            0, & \text{otherwise}
        \end{cases}

:math:`C_{i,j}` is the square adjacency matrix containing the distances (weights)
between city i and city j, with i and j :math:`\in = \{1, 2, .., n\}`. 
If the problem is symmetric, then the matrix is symmetric.
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

The dimensions of the constraint and of the chromosome grow quadratically with the number of cities. This makes the FULL encoding
particularly unefficient in PaGMO when large dimensions are involved.

TSP encoded using Random Keys (RANDOMKEYS encoding)
---------------------------------------------------

In this encoding the sequence of the cities visited is obtained from a sequence of floats by considering their indices after sorting.

.. math::

        [0.23,0.54,0.01,0.6] \rightarrow [2,0,1,3]

This maps the TSP into a box-bounded continuous problem getting rid of the constrains and of all the problems associated with
the possibility of generating unfeasible offsprings.

For more information see: Bean, J. C. (1994). Genetic algorithms and random keys for sequencing and optimization. ORSA journal on computing, 6(2), 154-160.

TSP encoded using the city sequence (CITIES encoding)
-----------------------------------------------------

In this encoding the sequence of the cities visited is directly put in the chromosome, creating a problem with one constraint: the
sequence must be a permutation of the city indexes.

.. math::

        [0,3,6,4,2,1,5] \rightarrow VALID
        
        [5,2,6,4,3,2,5] \rightarrow INVALID
