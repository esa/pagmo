.. _hypervolume:

Hypervolume
===========

A Quick Look
------------

Many multiple objective optimization problems often require a good method for comparison of solutions.
The question of how to do it is neither trivial, nor definite.
Hypervolume indicator (also known as Lebesgue measure or S-metric) found its application in that domain.
PyGMO allows for computing the hypervolume for the population objects (in that case, each individual fitness vector is treated as d-dimensional point), or for the 
fixed set of points.

Main class for that purpose is the :class:`PyGMO.util.hypervolume` class.
It is a main interface for computation of the hypervolume.
In order to compute a hypervolume from a population object, you can do the following:

.. code-block:: python

   from PyGMO import *
   from PyGMO.util import *
   prob = problem.dtlz(1)
   pop = population(prob,20)
   hv = hypervolume(pop)
   print hv.compute(r = [1000.0]*3)

Additionally, you can use the same module for the computation of a fixed hypervolume problem:

.. code-block:: python

   from PyGMO.util import *
   hv = hypervolume([[3,1],[2,1], [3,1]])
   print hv.compute(r = [4.0]*2)

You can also explicitly request an particular algorithm for the computation:

.. code-block:: python

   from PyGMO.util import *
   hv = hypervolume([[2,1,1] , [1,2,1], [1,1,2]])
   print hv.compute(r = [3.5]*3, algorithm = hv_algorithm.hv3d())

Besides the computation of the hypervolume of given pareto front (which may serve as a quality indicator), other methods are available as well:

.. code-block:: python

   from PyGMO.util import *
   hv = hypervolume([[2,1,1] , [1,2,1], [1,1,2]])
   print hv.exclusive(p_idx = 0, r = [3.5]*3)  # returns the exclusive volume by point 0
   print hv.least_contributor(r = [3.5]*3)  # returns the index of the least contributor

Method 'exlusive' computes the hypervolume of the point at given index.
Method 'least_contributors' returns the index of the individual (point) contributing the least to the hypervolume.


Available Hypervolume algorithms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
================================= ========================================== ================================================================
Common Name                       Name in PyGMO                              Comments
================================= ========================================== ================================================================
hv2d                              :class:`PyGMO.util.hv_algorithm.hv2d`      Algorithm for quick computation of 2d hypervolumes.
hv3d                              :class:`PyGMO.util.hv_algorithm.hv3d`      Algorithm for quick computation of 3d hypervolumes.
hv4d                              :class:`PyGMO.util.hv_algorithm.hv4d`      Algorithm for quick computation of 4d hypervolumes.
WFG                               :class:`PyGMO.util.hv_algorithm.wfg`       General algorithm for any dimension.
HOY                               :class:`PyGMO.util.hv_algorithm.hoy`       General algorithm for any dimension.
Bringmann-Friedrich approximation :class:`PyGMO.util.hv_algorithm.bf_approx` Algorithm for computation of the approximated least contributor.
FPRAS                             :class:`PyGMO.util.hv_algorithm.bf_fpras`  Algorithm for computation of the approximated hypervolume.
================================= ========================================== ================================================================

Detailed Documentation
----------------------
.. class:: PyGMO.util.hypervolume()

   This class allows for setting up hypervolume computation problems.
   Given hypervolume problem can be set up using population object, or by a list object.

   .. method:: __init__((PyGMO.population)pop)

      Constructs a hypervolume problem from a population object.
      In that case, each individual's fitness vector is pulled from the population, and treated as a point
      in hyperspace.

      USAGE:
         from PyGMO import *
         from PyGMO.util import *

         prob = problem.dtlz(1)

         pop = population(prob,20)

         hv = hypervolume(pop)

   .. method:: __init__((list)L)

      Constructs a custom hypervolume problem from a list.
      List object must contain other list objects that represent points in hyperspace.
      List object cannot be empty, and the dimension of each point must be no lesser than 2.

      USAGE:
         from PyGMO.util import *

         hv = hypervolume([[2,1,1], [1,1,2], [1,2,1]])

   .. method:: compute(r, algorithm = None)

      Computes the hypervolume for given problem, using the provided reference point r.
      Keyword `algorithm` must be an instance of algorithms that can be found inside `PyGMO.util.hv_algorithm` module.
      If the keyword is not provided, PyGMO chooses one automatically using the information about the reference point.
      In case of 2, 3 and 4 dimensions, algorithms hv2d, hv3d and hv4d are used.
      For larger dimensions the default method is the WFG.
      As of yet, it is required that reference point is numerically no lesser by each dimension than any point from the previously constructed set of points.

      * r - reference point used for computation
      * algorithm (optional) - hypervolume algorithm used for the computation, uses the best performing algorithm for given dimension by default

      USAGE:
         print hv.compute([3,3,3])

         print hv.compute([3,3,3], algorithm = hv_algorithm.hv3d())

         print hv.compute([3,3,3], algorithm = hv_algorithm.wfg())

   .. method:: exclusive(p_idx, r, algorithm = None)
      
      Computes the exlusive hypervolume for point at given index 'p_idx', using the provided reference point 'r' and the hypervolume algorithm (optional).
      Keyword `algorithm` must be an instance of algorithms that can be found inside `PyGMO.util.hv_algorithm` module.
      If the keyword is not provided, PyGMO chooses one automatically using the information about the reference point.

      * p_idx - index of the point for which we compute the exclusive hypervolume
      * r - reference point used for computation
      * algorithm (optional) - hypervolume algorithm used for the computation, uses the best performing algorithm for given dimension by default

      USAGE:
         hv.exclusive(p_idx=5, r=[5.0]*3)

         hv.exclusive(p_idx=5, r=[5.0]*3, algorithm=hv_algorithm.hv3d())

   .. method:: least_contributor(r, algorithm = None)
      
      Computes the least contributor to the hypervolume using provided reference point 'r' and the hypervolume algorithm (optional).
      Keyword `algorithm` must be an instance of algorithms that can be found inside `PyGMO.util.hv_algorithm` module.
      If the keyword is not provided, PyGMO chooses one automatically using the information about the reference point.

      * r - reference point used for computation
      * algorithm (optional) - hypervolume algorithm used for the computation, uses the best performing algorithm for given dimension by default

      USAGE:
         hv.least_contributor(r=[5.0]*3)

         hv.least_contributor(r=[5.0]*3, algorithm=hv_algorithm.bf_approx())

.. class:: PyGMO.util.hv_algorithm.hv2d()

    This is the quick algorithm the 2 dimensional problems.

   .. method:: __init__()

      Creates an instance of `PyGMO.util.hv_algorithm.hv2d` class that serves as a parameter to the hypervolume object.

.. class:: PyGMO.util.hv_algorithm.hv3d()

    This is the quick algorithm for the 3 dimensional problems.

   .. method:: __init__()

      Creates an instance of `PyGMO.util.hv_algorithm.hv3d` class that serves as a parameter to the hypervolume object.

.. class:: PyGMO.util.hv_algorithm.hv4d()

    This is the quick algorithm for the 4 dimensional problems.

   .. method:: __init__()

      Creates an instance of `PyGMO.util.hv_algorithm.hv4d` class that serves as a parameter to the hypervolume object.

.. class:: PyGMO.util.hv_algorithm.wfg()

    This is the implementation of the WFG algorithm.
    Its main purpose is handling hypervolume computation for any dimension.

   .. method:: __init__()

      Creates an instance of `PyGMO.util.hv_algorithm.wfg` class that serves as a parameter to the hypervolume object.

.. class:: PyGMO.util.hv_algorithm.hoy()

    This is the implementation of the HOY algorithm.
    Its main purpose is handling hypervolume computation for any dimension.

   .. method:: __init__()

      Creates an instance of `PyGMO.util.hv_algorithm.hoy` class that serves as a parameter to the hypervolume object.

.. class:: PyGMO.util.hv_algorithm.bf_approx()

    This is the implementation of the Bringmann-Friedrich approximation algorithm.
    Its main purpose is handling least contributor computation for any dimension.
    Algorithm's output is an approximation of the exact one to a certain degree of accuracy and certain confidence.

   .. method:: __init__(use_exact = True, trivial_subcase_size = 1, eps = 1e-1, delta = 1e-4, gamma = 0.25, delta_multiplier = 0.775, initial_delta_coeff = 1e-1, alpha = 0.2)

      Creates an instance of `PyGMO.util.hv_algorithm.bf_approx` class that serves as a parameter to the hypervolume object.
      Default values for the parameters of the algorithm were obtained from the shark implementation of the algorithm:
      http://image.diku.dk/shark/doxygen_pages/html/_least_contributor_approximator_8hpp_source.html

      Parameters:
      	* use_exact - should bf_approx use exact methods for computation
      	* trivial_subcase_size - when the number of points overlapping the bounding box is smaller or equal to that argument, we compute the exlusive hypervolume exactly
      	* eps - accuracy of approximation
      	* delta - confidence of approximation
      	* gamma - constant used for computation of delta for each of the points during the sampling
      	* delta_multiplier - factor with which delta diminishes each round
      	* initial_delta_coeff - initial coefficient multiplied by the delta at round 0
      	* alpha - coefficicient stating how accurately current lowest contributor should be sampled

.. class:: PyGMO.util.hv_algorithm.bf_fpras()

    This is the implementation of the Bringmann-Friedrich FPRAS algorithm, applied to the hypervolume computation problem.
    Its main purpose is handling hypervolume computation for any dimension.
    Algorithm's output is an approximation of the exact one to a certain degree of accuracy and certain confidence.

   .. method:: __init__(eps = 1e-2, delta = 1e-2)

      Creates an instance of `PyGMO.util.hv_algorithm.bf_fpras` class that serves as a parameter to the hypervolume object.

      Parameters:
      	* eps - accuracy of approximation
      	* delta - probability of error of the approximation
