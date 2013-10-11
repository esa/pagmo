.. _racing_algorithms:

=======================================================================
Racing algorithms
=======================================================================

Given an optimization problem, the solution provided by a stochastic optimizers
will not be deterministic, as their names implied. Nevertheless, the intrinsic
properties of each individual optimizers may make themselves perform better on
some classes of problems and worse on certain problems. Well, we are not too
concerned about this No Free Lunch Theorem here. There is another interesting
question: Given a problem (or a list of problems), which optimization
algorithms at our disposal is the most appropriate one, that will yield the
best optimization performance? Is there a way to find the answer without
resorting to pure brute force approach?

The approach of racing algorithms is based on statistical testing. Each
algorithm will be evaluated iteratively throughout the race. The performance of
an algorithm is defined as its ability to evolve a population, measured by the
fitness of the evolved champion. Whenever statistical significant evidence
about the inferiority of an algorithm is found, that algorithm will be dropped
out of the race. Naturally, it is guaranteed that the eventually survived
algorithms are statistically significantly better than those which are
discarded. The underlying statistical testing machinery used is Friedman test.

Using ``race_algo`` to race algorithms
######################################

Let us first set up a a problem and the algorithms of interest.

.. code-block:: python

    from PyGMO import *
    prob = problem.ackley()
    algo_list = []
    algo_list.append(algorithm.de(gen=100))
    algo_list.append(algorithm.cmaes(gen=100))
    algo_list.append(algorithm.pso(gen=100))

``race_algo`` contains the mechanisms that can be used to race different
algorithms. Its interface is very similar to that of ``pop_race``. The
following code sets up the necessary racing object and start the race. The
population size intialization argument specifies how large the internal
populations to be evolved should be, while each algorithm is being evaluated.

.. code-block:: python

    racer = util.race_algo(algo_list, prob, pop_size=100, seed=0)
    winners, n_evaluation = racer.run(n_final=1)

There are a couple of arguments that can be specified for ``run``, such as
``max_count`` which indicates the allowed computation budget. As above,
the indices of the winning algorithm are stored in ``winners``, while
``n_evaluation`` reflects how many calls to ``evolve()`` has been made.

In addition, it is possible to cover a range of problems simultaneously.

.. code-block:: python

    probs = [problem.ackley(), problem.griewank()]
    racer = util.race_algo(algo_list, probs, pop_size=100, seed=0)
    winners, n_evaluation = racer.run(n_final=1)

Internally, when ``race_algo`` is assigning a performance measure to each
algorithm, it randomly samples one of the problems and evolve popualtion with
respect to that problem.

Example: Identifying good configurations of an algorithm
########################################################

This example shows how ``race_algo`` can be used to select the good
configurations of an algorithm. Taking ``pso_gen`` as an example, there are
quite a few parameters to be configured. We will focus on the three of the
parameters as shown in the following code.

.. code-block:: python

    from PyGMO import *
    import itertools

    prob = problem.ackley()

    pop_size = 20
    fits = []

    variants = [1,2,3,4,5,6]
    neighb_types = [1,2,3,4]
    neighb_param = [1,2,3,4,5]
    pars = list(itertools.product(variants, neighb_types, neighb_param))
    algo_list = []
    for p in pars:
        args = {}
        args['variant'], args['neighb_type'], args['neighb_param'] = p
        algo_list.append(algorithm.pso(gen=100,**args))

    racer = util.race_algo(algo_list, prob, pop_size=pop_size, seed=0)
    winners, n_evaluated = racer.run(5, 0, len(algo_list)*30, 0.05, [], True, True)

    print 'Winning algorithm:'
    for i in winners:
        print zip(['variant', 'neighb_type', 'neighb_param'], pars[i])
        
    print 'Evaluated algorithms for %d times' % n_evaluated

    algo_default = algorithm.pso_gen(gen=100)
    algo_selected = algo_list[winners[0]]
    pop = population(prob, pop_size)
    pop_default = algo_default.evolve(pop)
    pop_selected = algo_selected.evolve(pop)
    print 'Default algorithm:', pop_default.champion.f
    print 'Selected algorithm:', pop_selected.champion.f

The output should be similar to the following (feel free to try a few more times due to the effects of random population initialization):

.. code-block:: python

    Winning algorithm:
    [('variant', 5), ('neighb_type', 1), ('neighb_param', 1)]
    [('variant', 5), ('neighb_type', 1), ('neighb_param', 2)]
    [('variant', 5), ('neighb_type', 1), ('neighb_param', 3)]
    [('variant', 5), ('neighb_type', 1), ('neighb_param', 4)]
    [('variant', 5), ('neighb_type', 1), ('neighb_param', 5)]
    Evaluated algorithms for 354 times
    Default algorithm: (0.8463040459318276,)
    Selected algorithm: (0.05200250235673254,)

Initially, the combinations of parameters result in a total of 120
configurations to be chosen from. Using ``race_pop``, staistical significant
results can be obtained after 354 evaluations. Without racing mechanism, this
is translated to running about only 3 trials on each algorithm configuration,
based on which statistical significant conclusion is difficult to be established.
This issue, as well as the appropriate distribution of computational budget, is
automatically handled by ``race_algo``.

Why not choosing a single best winner?
**************************************
However, ``race_algo`` is not magic. Consider changing ``n_final`` in the above
example from 5 to 4 (or less). The output are:

.. code-block:: python

    Winning algorithm:
    [('Variant', 5), ('neighb_type', 1), ('neighb_param', 1)]
    [('Variant', 5), ('neighb_type', 1), ('neighb_param', 2)]
    [('Variant', 5), ('neighb_type', 1), ('neighb_param', 3)]
    [('Variant', 5), ('neighb_type', 1), ('neighb_param', 4)]
    Evaluated algorithms for 3599 times
    Default algorithm: (2.187472414230023,)
    Selected algorithm: (0.04567759094164403,)

The consumed evaluation count is dramatically larger! The reason is that the
last 5 surviving algorithms are so close to each other that it is very
difficult to statistically significantly distinguish their performance. In fact, the race is forcefully terminated as the provided budget is exhausted.

.. note::
    Racing mechanism is most useful when the diversity of the entities in the
    racing pool is large, as the worse ones can be quickly weeded out with high
    staistical confidence, saving the evaluation budget the rest. When the
    entities are too similar to each other, the required evaluation will be
    very high in order to derive statistical significant conclusion.

