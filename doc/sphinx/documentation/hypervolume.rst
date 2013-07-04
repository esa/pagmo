Hypervolume
==========

A Quick Look
------------

Many multiple objective optimization problems often require a good method for comparison of solutions.
The question of how to do it is neither trivial, nor definite.
Hypervolume indicator (also known as Lebesgue measure or S-metric) found its application in that domain.
PyGMO allows for computing the hypervolume for the population objects (in that case, each individual's fitness vector is treated as d-dimensional point), or for the 
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
Lebmeasure                         :class:`PyGMO.hv_algorithm.lebmeasure`   General algorithm for any number of dimensions > 2, (very slow!)
================================== ======================================== ===================================================================

Detailed Documentation
----------------------
.. autoclass:: PyGMO.hypervolume._base()

   .. automethod:: PyGMO.hypervolume._base.compute

.. autoclass:: PyGMO.hypervolume.base()

.. autoclass:: PyGMO.hv_algorithm.lebmeasure
