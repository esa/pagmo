.. _tsp:

================================================================
Traveling Salesman Problems with PyGMO
================================================================

.. toctree::
   :maxdepth: 2

   create_tsp
   

=======================================================================
Traveling Salesman Problems
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

