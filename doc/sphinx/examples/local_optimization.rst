.. _adding_a_new_algorithm:

===============================
Nonlinear optimization in PyGMO
===============================

Besides a variety of available metaheuristics and global optimization algorithms, PaGMO interfaces several popular local optimizers (full list available on page :ref:`algorithms`).
Local optimization algorithms are interfaced such that they are compatible with all of other PyGMO features, such as Python interface, island model or algorithm evaluation metrics.
In this example we will present several features and usecases of nonlinear local optimizers available in PyGMO.

.. note::
 Almost all of the local optimizers (with an exception of :class:`PyGMO.algorithm.cs` available in PaGMO are interfaces to third-party dependencies.
 When using Python local optimizers (such as :class:`PyGMO.algorithm.scipy_fmin`), corresponding Python dependencies have to be installed.
 To use other high-performance optimizers (e.g.: :class:`PyGMO.algorithm.snopt`), PaGMO has to be compiled with the corresponding flag (e.g.: ``ENABLE_SNOPT``), while the library has to be installed and visible to the build system.

In the remainder of this section we will show how to use the SNOPT and WHORP nonlinear optimizers in PyGMO.

SNOPT library
=============

Installation
------------

To use SNOPT you first need to obtain either the full version of the library, or request a trial version from here:
https://ccom.ucsd.edu/~optimizers/downloads.php

Once SNOPT is installed in the system, you need to recompile PaGMO with ``ENABLE_SNOPT=ON`` flag.

Usage
-----

The following script optimizes the :class:`PyGMO.problem.luksan_vlcek_1` problem using default SNOPT configuration.

.. code-block:: python

    from PyGMO import *

    prob = problem.luksan_vlcek_1(10)
    # One 'individual', i.e. initial condition, is sufficient for local optimization
    pop = population(prob, 1)  # Population with one individual
    alg = algorithm.snopt()  # Use SNOPT library
    pop = alg.evolve(pop)
    print(pop.champion.x)
    print(pop.champion.f)

Additionally, several parameters of SNOPT are exposed to the Python interface:

.. code-block:: python

    # Set the upper limit on number of iterations, feasibility and optimality tolerances
    # and print SNOPT output the screen.
    alg = algorithm.snopt(major_iter=1000, feas_tol=1e-10, opt_tol=1e-10, screen_output=True)
    pop = population(prob, 1)
    pop = alg.evolve(pop)
    print(pop.champion.x)
    print(pop.champion.f)

WORHP library
=============

Installation
------------
.. note::
 PyGMO is supporting WORHP with version 1.8

Compiled WORHP library interface can be found at http://www.worhp.de/.
The software package contains the binaries and the header files, which have to be installed in the system (see installation guide available on the official WORHP website).

WORHP optimizer, besides requiring third party dependencies, also checks for a valid license and a configuration file.
License file ``worhp.lic`` can either reside in the current working directory or be referenced by the ``WORHP_LICENSE_FILE`` environment variable (either full or absolute path with the filename). Similarly the parameter configuration file ``param.xml`` can reside in the working directory or be referenced by ``WORHP_PARAM_FILE``. PaGMO needs to be compiled witht the flag ``ENABLE_WORHP=ON``.

Usage
-----

Usecase of :class:`PyGMO.algorithm.worhp` is the same as that of the SNOPT as show previously.
Additionally to three key parameters exposed in the constructor, it is also possible to set some WORHP parameters through the setter.

.. note::
 Three parameters are already exposed in the constructor of worhp algorithm (``TolFeas``, ``TolOpti`` and ``MaxIter``).
 Each of those is also set to the default values other than the specified in the default configuration file.
 **All parameters set by a setter method, as well as the parameters in the constructor (also the defaults!) are oveloading the values from the configuration file ``param.xml``.**
 List of setter-exposed WORHP parameters is available at: https://github.com/esa/pagmo/blob/master/src/algorithm/worhp_cpp_wrapper/set_get_params.cpp

Script below shows a usecase of WORHP as local optimizer in PaGMO:

.. code-block:: python

    from PyGMO import *

    prob = problem.py_pl2pl()
    pop = population(prob, 1)
    alg = algorithm.worhp(MaxIter=1000, TolFeas=1e-10)
    alg.set_param("AcceptTolOpti", 1e-4)
    pop = alg.evolve(pop)
    print(pop.champion.x)
    print(pop.champion.f)

SNOPT and WORHP comparison
==========================
This simple script seves as a comparison of two algorithms when run on several basic benchmarks

First, we will define some functions later used for the experiment:

.. code-block:: python

    from PyGMO import population, problem, algorithm
    import PyKEP
    import sys
    import random

    def run(seed, alg, prob):
        pop = population(prob, 1, seed)
        pop = alg.evolve(pop)
        feas = prob.feasibility_x(pop.champion.x)
        return (pop.champion.f[0], pop.problem.fevals, feas)

    def plot(l, title):
        import matplotlib.pyplot as plt
        a1, a2 = zip(*l)
        plt.xscale('log')

        d = filter(lambda x: x[2], a1)
        if d:
            plt.scatter(zip(*d)[0], zip(*d)[1], c='r', s=30, marker='o',
                label='WORHP, {}/{} feasible'.format(len(d), len(a1)))
        d = filter(lambda x: x[2], a2)
        if d:
            plt.scatter(zip(*d)[0], zip(*d)[1], c='g', s=30, marker='o',
                label='SNOPT, {}/{} feasible'.format(len(d), len(a1)))
        plt.legend()
        plt.xlabel('Objective function value')
        plt.grid(True)
        plt.ylabel('Function evaluations')
        plt.title(title)
        plt.show()

We compute 100 solutions to Luksan-Vlcek problem (dim=10) using SNOPT and WORHP, and plot the obtained objective function value against the number of function evaluations. For each comparison we start with the same initial conditions.
Results of the script below will open a new matpotlib window:

.. code-block:: python

    prob = problem.luksan_vlcek_1(10)

    l = []
    N = 100  # number of restarts
    sys.stdout.write("Computing {} luksan_vlcek_1 solutions..\n".format(N))

    for i in range(1, N + 1):
        tmp = []
        seed = random.randint(0,1000000)  # We start from the same initial conditions
        tmp.append(run(seed, algorithm.worhp(MaxIter=1000, TolFeas=1e-8, TolOpti=1e-8), prob))
        tmp.append(run(seed, algorithm.snopt(major_iter=1000, feas_tol=1e-8, opt_tol=1e-8), prob))
        l.append(tmp)
        sys.stdout.write("{}/{} ".format(i, N))
    plot(l, "Luksan Vlcek 1 (D=10)")


Similarly we compute 20 solutions to Python planet to planet problem (PyGMO needs to be compiled with ``ENABLE_GTOC=ON``):

.. code-block:: python

    # Warning: This is a very long test, Change N to something smaller (5) to get a quick overview
    N = 20  # number of restarts

    for n_seg in [5, 10, 20, 40]:
        l = []
        sys.stdout.write("Computing {} py_pl2pl solutions for {} segments..\n".format(N, n_seg))
        for i in range(1, N + 1):
            prob = problem.py_pl2pl(nseg=5)
            tmp = []    
            seed = random.randint(0,1000000)  # We start from the same initial conditions
            tmp.append(run(seed, algorithm.worhp(MaxIter=1000, TolFeas=1e-8, TolOpti=1e-8), prob))
            tmp.append(run(seed, algorithm.snopt(major_iter=1000, feas_tol=1e-8, opt_tol=1e-8), prob))
            l.append(tmp)
            sys.stdout.write("{}/{} ".format(i, N))
        plot(l, "Python Planet to Planet (N segments={})".format(n_seg))
