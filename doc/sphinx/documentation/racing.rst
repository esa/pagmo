Racing
======

Overview
--------
Two types of racing are supported in PyGMO, namely the racing between
individuals in a population, and racing between algorithms.

Racing of individuals in a population can be done as the following. Note that
in *PyGMO.population*, there is also a *race* method which can be conveniently
called. However, in order to utilize the caching mechanism in racing, one has
to invoke racing explicitly via *PyGMO.util.race_pop*.

.. code-block:: python

    from PyGMO import *
    prob = problem.noisy(problem.ackley(2))
    pop = population(prob,20)
    racer = util.race_pop(population=pop,seed=0)
    winners, fevals = racer.run(1) # Extract the single best winner

Illustration of the effects of caching:

.. code-block:: python

    from PyGMO import *
    pop = population(problem.noisy(), 10, 0)

    # Race overlapping individuals without caching
    winners, fevalsA = pop.race(1, racers_idx=range(0,5))
    winners, fevalsB = pop.race(1, racers_idx=range(3,8))

    # Race overlapping individuals with caching via race_pop
    racer = util.race_pop(pop)
    winners, fevalsA_cache = racer.run(1, racers_idx=range(0,5))
    winners, fevalsB_cache = racer.run(1, racers_idx=range(3,8))

    print 'Consumed fevals', fevalsA + fevalsB                          # 39
    print 'Consumed fevals [caching]', fevalsA_cache + fevalsB_cache    # 26

Racing of algorithms can be done as follows:

.. code-block:: python

    from PyGMO import *
    algos = [algorithm.de(), algorithm.cmaes(), algorithm.pso()]

    prob = problem.ackley()                            # To race over a single problem
    racer1 = util.race_algo(algos, prob)
    winners1, fevals1 = racer1.run(1)

    probs = [problem.ackley(), problem.griewank()]     # To race over multiple problems
    racer2 = util.race_algo(algos, probs)
    winners2, fevals2 = racer2.run(1)

Detailed documentation
----------------------

.. autoclass:: PyGMO.util.race_pop
    :members:
    
    .. automethod:: PyGMO.util.race_pop.__init__()

.. autoclass:: PyGMO.util.race_algo
   :members:

   Similar to race_pop, this class contains the mechanisms to race different
   algorithms. The race among the algorithms can be performed over a single
   problem, or over a multiple of them.

   .. automethod:: PyGMO.util.race_algo.__init__()
