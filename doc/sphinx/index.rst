========================================
Welcome to PyGMO (V1.0)!
========================================

PyGMO (the Python Parallel Global Multiobjective Optimizer) is a scientific library providing a large number
of optimisation problems and algorithms under the same powerful parallelization
abstraction built around our 'generalised island-model' paradigm. What this means to the user is that the available algorithms
are all automatically parallelized (asynchronously, coarse-grained approach) thus making efficient use of the underlying multicore
architecture. The user can also program his own solvers ... they also will be parallelized by PyGMO!! PyGMO's
implementation of the migration operator allows the user to easily define "migration paths" between a large number of "islands" (CPU cores).

Efficient implementantions of state-of-the-art bio-inspired algorithms are sided to state-of the art optimization algorithms (Simplex Methods, SQP methods ....)
and can be easily mixed (also with your newly invented algorithms) to build a super-algorithm exploiting cooperation via the asynchronous, generalized island model.

PyGMO can be used to solve constrained, unconstrained, single objective, multiple objective, contiuous, mixed int 
optimization problem, or to perform research on novel algorithms and paradigms and easily compare them to state of the art
implementations of established ones.

PyGMO is interfaced with SciPy optimization algorithms, NLOPT algorithms, SNOPT, IPOPT and, hopefully .... more to come.

NOTE: This is a pre-release web-site for PyGMO. Please send your comments via the `pagmo/PyGMO mailing list <http://sourceforge.net/mail/?group_id=238743>`_ and submit any
bugs via our `pagmo/PyGMO bug-tracker <http://sourceforge.net/tracker/?group_id=238743&atid=1133009>`_

Install
=======

You may download the latest code from our `git repo <http://sourceforge.net/scm/?type=git&group_id=238743>`_
To compile ... select, using ccmake, the option compile PyGMO, then make and install

.. toctree::
   :maxdepth: 2

   pre_requisites
   download
   installation

Quick-Start
===========

.. toctree::
   :maxdepth: 2

   quickstart

Documentation
=============

.. toctree::
   :maxdepth: 2

   documentation/documentation

Examples
=============

.. toctree::
   :maxdepth: 2

   examples

Credits
=======

.. toctree::
   :maxdepth: 2

   credits
   

