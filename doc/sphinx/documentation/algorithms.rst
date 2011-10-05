Algorithms
==========

A Quick Look
------------

Algorithms in PyGMO are constructed and then used to evolve a population via their evolve method. The user can
implement his own algorithm in Python (in which case they need to derive from :class:`PyGMO.algorithm.base`)

================================== ========================================= =============== ============== ================ ===========================================
Common Name                        Name in PyGMO                             Chromosome      Constraints    Fitness          Comments
================================== ========================================= =============== ============== ================ ===========================================
Differential Evolution (DE)        :class:`PyGMO.algorithm.de`               Continuous      Unconstrained  Single-objective 
Particle Swarm Optimization (PSO)  :class:`PyGMO.algorithm.pso`              Continuous      Unconstrained  Single-objective 
Particle Swarm Optimization (PSO)  :class:`PyGMO.algorithm.pso_gen`          Continuous      Unconstrained  Single-objective also problems deriving from base_stochastic
Simple Genetic Algorithm (SGA)     :class:`PyGMO.algorithm.sga`              Continuous      Unconstrained  Single-objective 
Corana's Simulated Annealing (SA)  :class:`PyGMO.algorithm.sa_corana`        Continuous      Unconstrained  Single-objective 
Artificial Bee Colony (ABC)        :class:`PyGMO.algorithm.bee_colony`       Continuous      Unconstrained  Single-objective 
Firefly Optimization (FF)          :class:`PyGMO.algorithm.firefly`          Continuous      Unconstrained  Single-objective 
Monotonic Basin Hopping (MBH)      :class:`PyGMO.algorithm.mbh`                    N/A            N/A              N/A       
Multistart (MS)                    :class:`PyGMO.algorithm.ms`                     N/A            N/A              N/A       
Improved Harmony Search (IHS)      :class:`PyGMO.algorithm.ihs`              Mixed-Int       Unconstrained  Single-objective 
Compass Search (CS)                :class:`PyGMO.algorithm.cs`               Continuous      Unconstrained  Single-objective  
Monte Carlo (MC)                   :class:`PyGMO.algorithm.monte_carlo`      Mixed-Int       Constrained    Single-objective  
Monte Carlo (MC)                   :class:`PyGMO.algorithm.py_example`       Mixed-Int       Constrained    Single-objective 
Ant Colony optimization (ACO)      :class:`PyGMO.algorithm.aco`                    Int       Unconstrained  Single-objective only problems deriving from base_aco
Cross-Entropy Method (CE)          :class:`PyGMO.algorithm.py_cross_entropy` Continuous      Unconstrained  Single-objective 
Nelder-Mead simplex                :class:`PyGMO.algorithm.scipy_fmin`       Continuous      Unconstrained  Single-objective requires scipy module installed
L-BFGS-B                           :class:`PyGMO.algorithm.scipy_l_bfgs_b`   Continuous      Unconstrained  Single-objective requires scipy module installed
Sequential Least SQuares Prog.     :class:`PyGMO.algorithm.scipy_slsqp`      Continuous      Constrained    Single-objective requires scipy module installed
================================== ========================================= =============== ============== ================ ===========================================

Detailed Documentation
----------------------
.. autoclass:: PyGMO.algorithm._base()

   .. automethod:: PyGMO.algorithm._base.evolve

.. autoclass:: PyGMO.algorithm.base()

.. autoclass:: PyGMO.algorithm.de

   .. automethod:: PyGMO.algorithm.de.__init__

.. autoclass:: PyGMO.algorithm.pso

   .. automethod:: PyGMO.algorithm.pso.__init__

.. autoclass:: PyGMO.algorithm.pso_gen

   .. automethod:: PyGMO.algorithm.pso_gen.__init__

.. autoclass:: PyGMO.algorithm.sga

   .. automethod:: PyGMO.algorithm.sga.__init__

   .. attribute:: PyGMO.algorithm.sga.mutation.RANDOM

     Random mutation (width is set by the width argument in :class:`PyGMO.algorithm.sga`)

   .. attribute:: PyGMO.algorithm.sga.mutation.GAUSSIAN

     Gaussian mutation (bell shape standard deviation is set by the width argument in :class:`PyGMO.algorithm.sga` multiplied by the box-bounds width)

   .. attribute:: PyGMO.algorithm.sga.selection.ROULETTE

     Roulette selection mechanism

   .. attribute:: PyGMO.algorithm.sga.selection.BEST20

     Best 20% individuals are inserted over and over again

   .. attribute:: PyGMO.algorithm.sga.crossover.BINOMIAL

     Binomial crossover

   .. attribute:: PyGMO.algorithm.sga.crossover.EXPONENTIAL

     Exponential crossover

.. autoclass:: PyGMO.algorithm.sa_corana

   .. automethod:: PyGMO.algorithm.sa_corana.__init__

.. autoclass:: PyGMO.algorithm.bee_colony

   .. automethod:: PyGMO.algorithm.bee_colony.__init__

.. autoclass:: PyGMO.algorithm.firefly(*args)

   .. automethod:: PyGMO.algorithm.firefly.__init__

.. autoclass:: PyGMO.algorithm.ms

   .. automethod:: PyGMO.algorithm.ms.__init__

   .. method:: PyGMO.algorithm.ms.screen_output((bool) active)

      Activates screen output. Suggested only when no parallel optimization are performed.

   .. attribute:: PyGMO.algorithm.ms.algorithm

      Algorithm to be multistarted

.. autoclass:: PyGMO.algorithm.mbh

   .. automethod:: PyGMO.algorithm.mbh.__init__

   .. method:: PyGMO.algorithm.mbh.screen_output((bool) active)

      Activates screen output. Suggested only when no parallel optimization are performed.

   .. attribute:: PyGMO.algorithm.mbh.algorithm

      Algorithm to perform mbh 'local' search

.. autoclass:: PyGMO.algorithm.cs

   .. automethod:: PyGMO.algorithm.cs.__init__

.. autoclass:: PyGMO.algorithm.ihs

   .. automethod:: PyGMO.algorithm.ihs.__init__

.. autoclass:: PyGMO.algorithm.monte_carlo

   .. automethod:: PyGMO.algorithm.monte_carlo.__init__

.. autoclass:: PyGMO.algorithm.py_example

   .. automethod:: PyGMO.algorithm.py_example.__init__

.. autoclass:: PyGMO.algorithm.py_cross_entropy

   .. automethod:: PyGMO.algorithm.py_cross_entropy.__init__

.. autoclass:: PyGMO.algorithm.scipy_fmin

   .. automethod:: PyGMO.algorithm.scipy_fmin.__init__

.. autoclass:: PyGMO.algorithm.scipy_l_bfgs_b

   .. automethod:: PyGMO.algorithm.scipy_l_bfgs_b.__init__

.. autoclass:: PyGMO.algorithm.scipy_slsqp

   .. automethod:: PyGMO.algorithm.scipy_slsqp.__init__
