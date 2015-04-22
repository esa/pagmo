.. _adding_a_new_algorithm:

================================================================
Local optimization in PyGMO
================================================================

Besides a variety of available metaheuristics and global optimization algorithms, PaGMO interfaces several popular local optimizers (full list available on page :ref:`algorithms`).
Local optimization algorithms are interfaced such that they are compatible with all of other PyGMO features, such as Python interface, island model or algorithm evaluation metrics.

.. note::
 Almost all of the local optimizers available in PaGMO are interfaces to third-party dependencies.
 When using Python local optimizers (such as :ref:`PyGMO.algorithm.scipy_fmin`), corresponding Python dependencies have to be installed.
 To use other high-performance optimizers (e.g.: :ref:`PyGMO.algorithm.snopt`), PaGMO has to be compiled with the corresponding flag (e.g.: ``ENABLE_SNOPT``), while the library has to be installed and visible to the build system.

In the remainder of this section we will show how to use the SNOTP and WHORP optimizers in PyGMO.

SNOPT library
=============

The following script optimizes the ``luksan_vlcek_1`` problem using default SNOPT configuration.

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

WORHP optimizer, besides requiring the third party dependencies, also requires a valid license and a configuration file.

.. note::
 License file ``worhp.lic`` can either reside in the current working directory or be referenced by the ``WORHP_LICENSE_FILE`` environment variable (either full or absolute path, along with the filename). Similarly the parameter configuration file ``param.xml`` can reside in the working directory or be referenced by ``WORHP_PARAM_FILE``.

Usecase of ``PyGMO.algorithm.worhp`` is the same as that of the SNOPT as show previously.
Additionally to several parameters exposed in the constructor, it is also possible to set some WORHP parameters through the setter.

.. note::
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
