Hypervolume
==========

A Quick Look
------------

Many multiple objective optimization problems often require a good method for comparison of solutions.
The question of how to do it is neither trivial, nor definite.
Hypervolume indicator (also known as Lebesgue measure or S-metric) found its application in that domain.
PyGMO allows for computing the hypervolume for the population objects (in that case, each individual fitness vector is treated as d-dimensional point), or for the 
fixed set of points.

Main class for that purpose is the :class:`PyGMO.hypervolume` class.
It is a main interface for computation of the hypervolume.
In order to compute a hypervolume from a population object, you can do the following:

.. code-block:: python

   from PyGMO import *
   prob = problem.dtlz1(fdim=3)
   pop = population(prob,20)
   hv = hypervolume(pop)
   print hv.compute(r = [1000.0]*3)

Additionally, you can use the same module for the computation of a fixed hypervolume problem:

.. code-block:: python

   from PyGMO import *
   hv = hypervolume([[3,1],[2,1], [3,1]])
   print hv.compute(r = [4.0]*2)

You can also explicitly request an particular algorithm for the computation:

.. code-block:: python

   from PyGMO import *
   hv = hypervolume([[2,1,1] , [1,2,1], [1,1,2]])
   print hv.compute(r = [3.5]*3, algorithm = hv_algorithm.optimal3d())

Available Hypervolume algorithms
^^^^^^^^^^^^^^^^^^^^^^
================================== ======================================== ===================================================================
Common Name                        Name in PyGMO                            Comments
================================== ======================================== ===================================================================
Optimal2D                          :class:`PyGMO.hv_algorithm.optimal2d`    Algorithm for optimal computation of 2D hypervolumes.
Optimal3D                          :class:`PyGMO.hv_algorithm.optimal3d`    Algorithm for optimal computation of 3D hypervolumes.
LebMeasure                         :class:`PyGMO.hv_algorithm.lebmeasure`   General algorithm for any size of dimension > 2, (caution: very slow)
================================== ======================================== ===================================================================

Detailed Documentation
----------------------
.. class:: PyGMO.hypervolume()

   This class allows for setting up hypervolume computation problems.
   Given hypervolume problem can be set up using population object, or by a list object.

   .. method:: __init__((PyGMO.population)pop)

      Constructs a hypervolume problem from a population object.
      In that case, each individual's fitness vector is pulled from the population, and treated as a point
      in hyperspace.

      USAGE:
         from PyGMO import *

         prob = problem.dtlz1(fdim=3)

         pop = population(prob,20)

         hv = hypervolume(pop)

   .. method:: __init__((list)L)

      Constructs a custom hypervolume problem from a list.
      List object must contain other list objects that represent points in hyperspace.
      List object cannot be empty, and the dimension of each point must be no lesser than 2.

      USAGE:
         from PyGMO import *

         hv = hypervolume([[2,1,1], [1,1,2], [1,2,1]])

   .. method:: compute(r, algorithm = None)

      Computes the hypervolume for given problem, using the provided reference point r.
      Keyword `algorithm` must be an instance of algorithms that can be found inside `PyGMO.hv_algorithm` module.
      If the keyword is not provided, PyGMO chooses one automatically using the information about the reference point.
      In case of 2 and 3 dimensions, methods Optimal2D and Optimal3D are used.
      For larger dimensions the default method is the LebMeasure.
      As of yet, it is required that reference point is numerically no lesser by each dimension than any point from the previously constructed set of points.

      USAGE:
         print hv.compute([3,3,3])

         print hv.compute([3,3,3], algorithm = hv_algorithm.optimal3d())

         print hv.compute([3,3,3], algorithm = hv_algorithm.lebmeasure())

.. class:: PyGMO.hv_algorithm.optimal2d()

    This is the optimal algorithm the 2 dimensional problems.

   .. method:: __init__()

      Creates an instance of `PyGMO.hv_algorithm.optimal2d` class that serves as a parameter to the hypervolume object.

.. class:: PyGMO.hv_algorithm.optimal3d()

    This is the optimal algorithm for the 3 dimensional problems.

   .. method:: __init__()

      Creates an instance of `PyGMO.hv_algorithm.optimal3d` class that serves as a parameter to the hypervolume object.

.. class:: PyGMO.hv_algorithm.lebmeasure()

    This is the implementation of the LebMeasure algorithm.
    Its main purpose is handling hypervolume computation for 3 or more dimensions.
    This algorithm is quite slow due to its high computational complexity.

   .. method:: __init__()

      Creates an instance of `PyGMO.hv_algorithm.lebmeasure` class that serves as a parameter to the hypervolume object.
