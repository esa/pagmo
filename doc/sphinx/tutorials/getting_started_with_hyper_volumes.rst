.. _getting_started_with_hyper_volumes:

================================================================
Getting started with hypervolumes
================================================================

This tutorial will cover the features introduces by the hypervolumes.
First we will describe how the user interface to hypervolume computation was designed, and point out important notions that ought to be taken into account.
Later, we will go over several examples, in order to get you started with the hypervolume computation.

Hypervolume Interface and Construction
======================================

Main class used for the computation of the hypervolume indicator (also known as Lebesgue Measure or S-Metric) and other measures that derive from it, is the `PyGMO.util.hypervolume` class. Even though it is under the *util* package, importing everything from `PyGMO` using asterix, makes the hypervolume class accessible as if it would be right in `PyGMO` itself, i.e.:

.. code-block:: python

  from PyGMO.util import hypervolume
  'hypervolume' in dir()  # Returns True
    
.. code-block:: python

  from PyGMO import *
  'hypervolume' in dir()  # Returns True
  
Latter will the preferred method of importing the hypervolume interface. Since whole feature of the computation of hypervolume is bound tightly to many optimization algorithms, we provide two constructors for the hypervolume object:

.. code-block:: python

  from PyGMO import *

  prob = problem.dtlz2(fdim=3)  # Construct DTLZ-2 problem with 3-dimensional fitness space
  pop = population(prob, 50)  # Construct the population object
  hv = hypervolume(pop)  # Construct the hypervolume object from the population object
  
Since the fitness space of optimization problems have always been tied to the domain of computational geometry, first obvious constructor was to initiate a hypervolume problem out of a population object.
Inside the constructor, each point's fitness vector is copied, and added to the internal representation of the hypervolume problem.

In order to provide the users with a more general purpose engine for the computation of hypervolumes, it is also possible to supplement the constructor with a list or a tuple of points:

.. code-block:: python

  from PyGMO import *

  hv = hypervolume([[1,0],[0.5,0.5],[0,1]])

This constructor will initiate the hypervolume object with three, 2-dimensional points as the input point set. This type of explicit constructor is especially useful when the problem at hand is outside of the scope of the optimization framework, or when the direct constructor from population is too general (e.g. when we require some information on dominated points within given population).

Computing the hypervolume indicator and other features
======================================================

Before we give an overview of each hypervolume features, let us discuss the assumptions we make regarding the reference point, and the set of points in general:

1. We assume **minimization**, that is, a reference point is required to be numerically larger or equal in each objective, and strictly larger in at least one of them.
2. Although the hypervolume for one dimension is well defined mathematically, we require that any input data be at least 2-dimensional, be in input points or the reference point. Additionally, dimensions among the points and the dimension of reference point must be equal.


For simplicity, we will support every example below based on this simple 2-dimensional front:

.. code-block:: python

  from PyGMO import *

  hv = hypervolume( ((1, 0), (0.5, 0.5), (0, 1), (1.5, 0.75)) )
  r = (2,2)
  hv.compute(r)  # Returns 3.25 as an answer

We will refer to each point by it's position on X axis, e.g. first point is the point (0,1), fourth point is (1.5, 0.75) etc.

.. image:: ../images/tutorials/hv_front_2d_simple.png
  :width: 750px

Once the hypervolume object is created (either from the population object or a raw set of points), there are several measures we can request for:

1. **compute** - Returns the joint hypervolume of the set of points (S-Metric).

.. code-block:: python

  # hv and r refer to the data above
  hv.compute(r)  # Returns 3.25 as an answer

2. **exclusive** - Request for the computation of the exclusive hypervolume by point at given index (starting at 0).

.. code-block:: python

  # hv and r refer to the data above
  hv.exclusive(1, r)  # Returns 0.25 as an answer
  hv.exclusive(3, r)  # Returns 0.0 as an answer since third point is dominated

3. **least_contributor** - Returns the index of the point contributing the least to the hypervolume.

.. code-block:: python

  # hv and r refer to the data above
  hv.least_contributor(r)  # Returns 3 as an answer, since third point contributes no hypervolume

4. **greatest_contributor** - Returns the index of the point contributing the most to the hypervolume.

.. code-block:: python

  # hv and r refer to the data above
  hv.greatest_contributor(r)  # Returns either 0 or 2 as an answer

As you might have noticed, there is more than once candidate for the greatest contributor.
In cases like these, it is undefined which point will be returned. Similar non-policy applies for the least contributor (In practical cases, this scenario occurs even more often).

5. **contributions** - Returns a list of contributions for each of the point in a set. This may seem like a feature overloading the **exclusive** method, yet due to implementation details, explicit request for all contributions may be much faster (in the best case by the *O(n)* time).

.. code-block:: python

  # hv and r refer to the data above
  hv.contributions(r)  # Returns a tuple (0.5, 0.25, 0.5, 0.0)
