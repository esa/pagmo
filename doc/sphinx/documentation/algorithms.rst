.. _algorithms:

Algorithms
==========

A Quick Look
------------

Algorithms in PyGMO are objects, constructed and then used to optimize a problem via their evolve method. 
The user can implement his own algorithm in Python (in which case they need to derive from
:class:`PyGMO.algorithm.base`). You may follow the :ref:`adding_a_new_algorithm` tutorial. 
We also provide a number of algorithms that are considered useful for general purposes. 
Each algorithm can be associated only to problems of certain types: 
(Continuous, Integer or Mixed Integer)-(Constrained, Unconstrained)-(Single, Multi-objective).


Heuristic Optimization
^^^^^^^^^^^^^^^^^^^^^^
========================================= ========================================= =============== ===================================================================
Common Name                               Name in PyGMO                             Type            Comments
========================================= ========================================= =============== ===================================================================
Differential Evolution (DE)               :class:`PyGMO.algorithm.de`                    C-U-S      The original algorithm
Self-adaptive DE (jDE)                    :class:`PyGMO.algorithm.jde`                   C-U-S      Self-adaptive F, CR
DE with p-best crossover (mde_pbx)        :class:`PyGMO.algorithm.mde_pbx`               C-U-S      Self-adaptive F, CR
Differential Evolution (DE)               :class:`PyGMO.algorithm.de_1220`               C-U-S      Our own brew. self adaptive F, CR and variants 
Particle Swarm Optimization (PSO)         :class:`PyGMO.algorithm.pso`                   C-U-S      The PSO algorithm (canonical, with constriction factor, FIPS, etc.)
Particle Swarm Optimization (PSO)         :class:`PyGMO.algorithm.pso_gen`               C-U-S      Generational (also problems deriving from base_stochastic)
Simple Genetic Algorithm GRAY (SGA_GRAY)  :class:`PyGMO.algorithm.sga_gray`              C-U-S      Simple genetic algorithm with gray binary encoding
Simple Genetic Algorithm (SGA)            :class:`PyGMO.algorithm.sga`                  MI-U-S 
Vector Evaluated Genetic Algorithm (VEGA) :class:`PyGMO.algorithm.vega`                 MI-U-M      VEGA algorithm, multi-objective extension of SGA
(N+1)-EA Evol. Algorithm (SEA)            :class:`PyGMO.algorithm.sea`                   I-U-M      The multiobjective extension uses crowding distance operator
Non-dominated Sorting GA (NSGA2)          :class:`PyGMO.algorithm.nsga_II`               C-U-M      NSGA-II
S-Metric Selection EMOA (SMS-EMOA)        :class:`PyGMO.algorithm.sms_emoa`              C-U-M      Relies on the hypervolume computation.
Corana's Simulated Annealing (SA)         :class:`PyGMO.algorithm.sa_corana`             C-U-S 
Parallel Decomposition (PADE)             :class:`PyGMO.algorithm.pade`                  C-U-M      Parallel Decomposition (based on the MOEA/D framework)
Non-dominated Sorting PSO (NSPSO)         :class:`PyGMO.algorithm.nspso`                 C-U-M      Multi-Objective PSO
Strength Pareto EA 2 (SPEA2)              :class:`PyGMO.algorithm.spea2`                 C-U-M      Strength Pareto Evolutionary Algorithm 2
Artificial Bee Colony (ABC)               :class:`PyGMO.algorithm.bee_colony`            C-U-S 
Improved Harmony Search (IHS)             :class:`PyGMO.algorithm.ihs`                  MI-U-M      Integer and Multiobjetive not tested yet
Monte Carlo Search (MC)                   :class:`PyGMO.algorithm.monte_carlo`          MI-C-S
Monte Carlo Search (MC)                   :class:`PyGMO.algorithm.py_example`           MI-C-S      Written directly in Python
Covariance Matrix Adaptation-ES           :class:`PyGMO.algorithm.py_cmaes`              C-U-S      Written directly in Python
Covariance Matrix Adaptation-ES           :class:`PyGMO.algorithm.cmaes`                 C-U-S
========================================= ========================================= =============== ===================================================================

Meta-algorithms 
^^^^^^^^^^^^^^^
================================== ============================================ =============== ===========================================
Common Name                        Name in PyGMO                                Type            Comments
================================== ============================================ =============== ===========================================
Monotonic Basin Hopping (MBH)      :class:`PyGMO.algorithm.mbh`                    N/A          
Multistart (MS)                    :class:`PyGMO.algorithm.ms`                     N/A      
Augmented Lagrangian (AL)          :class:`PyGMO.algorithm.nlopt_auglag`          C-C-S         Requires PyGMO to be compiled with nlopt option. Minimization assumed
Augmented Lagrangian (AL)          :class:`PyGMO.algorithm.nlopt_auglag_eq`       C-C-S         Requires PyGMO to be compiled with nlopt option. Minimization assumed
Cstrs co-evolution                 :class:`PyGMO.algorithm.cstrs_co_evolution`    C-C-S         Minimization assumed
Cstrs Self-Adaptive                :class:`PyGMO.algorithm.cstrs_self_adaptive`   C-C-S         Minimization assumed
Cstrs Immune System                :class:`PyGMO.algorithm.cstrs_immune_system`   C-C-S         Immune system constraints handling technique
Cstrs CORE                         :class:`PyGMO.algorithm.cstrs_core`            C-C-S         CORE constraints handling technique (repairing technique)
================================== ============================================ =============== ===========================================

Local optimization 
^^^^^^^^^^^^^^^^^^
================================== ========================================= =============== =====================================================================
Common Name                        Name in PyGMO                             Type            Comments
================================== ========================================= =============== =====================================================================
Compass Search (CS)                :class:`PyGMO.algorithm.cs`                    C-U-S 
Nelder-Mead simplex                :class:`PyGMO.algorithm.scipy_fmin`            C-U-S      SciPy required. Minimization assumed
Nelder-Mead simplex                :class:`PyGMO.algorithm.gsl_nm`                C-U-S      Requires PyGMO to be compiled with GSL option. Minimization assumed
Nelder-Mead simplex variant 2      :class:`PyGMO.algorithm.gsl_nm2`               C-U-S      Requires PyGMO to be compiled with GSL option. Minimization assumed
Nelder-Mead simplex variant 2r     :class:`PyGMO.algorithm.gsl_nm2rand`           C-U-S      Requires PyGMO to be compiled with GSL option. Minimization assumed
Subplex (a Nelder-Mead variant)    :class:`PyGMO.algorithm.nlopt_sbplx`           C-C-S      Requires PyGMO to be compiled with nlopt option. Minimization assumed
L-BFGS-B                           :class:`PyGMO.algorithm.scipy_l_bfgs_b`        C-U-S      SciPy required. Minimization assumed
BFGS                               :class:`PyGMO.algorithm.gsl_bfgs2`             C-U-S      Requires PyGMO to be compiled with GSL option. Minimization assumed
BFGS 2                             :class:`PyGMO.algorithm.gsl_bfgs`              C-U-S      Requires PyGMO to be compiled with GSL option. Minimization assumed
Sequential Least SQuares Prog.     :class:`PyGMO.algorithm.scipy_slsqp`           C-C-S      SciPy required. Minimization assumed
Sequential Least SQuares Prog.     :class:`PyGMO.algorithm.nlopt_slsqp`           C-C-S      Requires PyGMO to be compiled with nlopt option. Minimization assumed
Truncated Newton Method            :class:`PyGMO.algorithm.scipy_tnc`             C-U-S      SciPy required. Minimization assumed
Conjugate Gradient (fr)            :class:`PyGMO.algorithm.gsl_fr`                C-U-S      Requires PyGMO to be compiled with GSL option. Minimization assumed
Conjugate Gradient (pr)            :class:`PyGMO.algorithm.gsl_pr`                C-U-S      Requires PyGMO to be compiled with GSL option. Minimization assumed
COBYLA                             :class:`PyGMO.algorithm.scipy_cobyla`          C-C-S      SciPy required. Minimization assumed
COBYLA                             :class:`PyGMO.algorithm.nlopt_cobyla`          C-C-S      Requires PyGMO to be compiled with nlopt option. Minimization assumed
BOBYQA                             :class:`PyGMO.algorithm.nlopt_bobyqa`          C-C-S      Requires PyGMO to be compiled with nlopt option. Minimization assumed
Method of Moving Asymptotes        :class:`PyGMO.algorithm.nlopt_mma`             C-C-S      Requires PyGMO to be compiled with nlopt option. Minimization assumed
SNOPT                              :class:`PyGMO.algorithm.snopt`                 C-C-S      Requires PyGMO to be compiled with snopt option. Minimization assumed
IPOPT                              :class:`PyGMO.algorithm.ipopt`                 C-C-S      Requires PyGMO to be compiled with ipopt option. Minimization assumed
================================== ========================================= =============== =====================================================================

Detailed Documentation
----------------------
.. class:: PyGMO.algorithm.base()

   All PyGMO algorithms derive from this class

   .. automethod:: PyGMO.algorithm.base.evolve

---------------

.. class:: PyGMO.algorithm.de

   .. automethod:: PyGMO.algorithm.de.__init__

---------------

.. class:: PyGMO.algorithm.jde

   .. automethod:: PyGMO.algorithm.jde.__init__

---------------

.. class:: PyGMO.algorithm.mde_pbx

   .. automethod:: PyGMO.algorithm.mde_pbx.__init__

---------------

.. class:: PyGMO.algorithm.de_1220

   .. automethod:: PyGMO.algorithm.de_1220.__init__

---------------

.. class:: PyGMO.algorithm.pso

   .. automethod:: PyGMO.algorithm.pso.__init__

---------------

.. class:: PyGMO.algorithm.pso_gen

   .. automethod:: PyGMO.algorithm.pso_gen.__init__

---------------

.. class:: PyGMO.algorithm.sea

   .. automethod:: PyGMO.algorithm.sea.__init__

---------------

.. class:: PyGMO.algorithm.sga

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

---------------

.. class:: PyGMO.algorithm.vega

   .. automethod:: PyGMO.algorithm.vega.__init__

   .. attribute:: PyGMO.algorithm.vega.mutation.RANDOM

     Random mutation (width is set by the width argument in :class:`PyGMO.algorithm.vega`)

   .. attribute:: PyGMO.algorithm.vega.mutation.GAUSSIAN

     Gaussian mutation (bell shape standard deviation is set by the width argument in :class:`PyGMO.algorithm.vega` multiplied by the box-bounds width)

   .. attribute:: PyGMO.algorithm.vega.crossover.BINOMIAL

     Binomial crossover

   .. attribute:: PyGMO.algorithm.vega.crossover.EXPONENTIAL

     Exponential crossover

---------------

.. class:: PyGMO.algorithm.sga_gray

   .. automethod:: PyGMO.algorithm.sga_gray.__init__

   .. attribute:: PyGMO.algorithm.sga_gray.mutation.UNIFORM

     Uniform mutation 

   .. attribute:: PyGMO.algorithm.sga_gray.selection.ROULETTE

     Roulette selection mechanism

   .. attribute:: PyGMO.algorithm.sga_gray.selection.BEST20

     Best 20% individuals are inserted over and over again

   .. attribute:: PyGMO.algorithm.sga_gray.crossover.SINGLE_POINT

     Single point crossover

---------------

.. class:: PyGMO.algorithm.nsga_II

   .. automethod:: PyGMO.algorithm.nsga_II.__init__

---------------

.. class:: PyGMO.algorithm.sms_emoa

   .. automethod:: PyGMO.algorithm.sms_emoa.__init__

---------------

.. class:: PyGMO.algorithm.pade

   .. automethod:: PyGMO.algorithm.pade.__init__
   
   .. attribute:: PyGMO.algorithm.pade.RANDOM

     Random generation of the weight vector
   
   .. attribute:: PyGMO.algorithm.pade.GRID

     Weight vectors are generated to equally divide the search space (requires a particular population size)

   .. attribute:: PyGMO.algorithm.pade.LOW_DISCREPANCY

     Weight vector are generated using a low discrepancy sequence

---------------

.. class:: PyGMO.algorithm.nspso

   .. automethod:: PyGMO.algorithm.nspso.__init__
   
   .. attribute:: PyGMO.algorithm.nspso.CROWDING_DISTANCE

     Individual with better crowding distance are prefered

   .. attribute:: PyGMO.algorithm.nspso.NICHE_COUNT

     Individuals with better niche count are prefered

   .. attribute:: PyGMO.algorithm.nspso.MAXMIN

     The MaxMin method is used to obtain the non-dominated set and to mantain diversity

---------------

.. class:: PyGMO.algorithm.spea2

   .. automethod:: PyGMO.algorithm.spea2.__init__

---------------

.. class:: PyGMO.algorithm.sa_corana

   .. automethod:: PyGMO.algorithm.sa_corana.__init__

---------------

.. class:: PyGMO.algorithm.bee_colony

   .. automethod:: PyGMO.algorithm.bee_colony.__init__

---------------

.. class:: PyGMO.algorithm.ms

   .. automethod:: PyGMO.algorithm.ms.__init__

   .. attribute:: PyGMO.algorithm.ms.screen_output

      When True, the algorithms produces output on screen 

   .. attribute:: PyGMO.algorithm.ms.algorithm

      Algorithm to be multistarted

---------------

.. class:: PyGMO.algorithm.mbh

   .. automethod:: PyGMO.algorithm.mbh.__init__

   .. attribute:: PyGMO.algorithm.mbh.screen_output

      When True, the algorithms produces output on screen 

   .. attribute:: PyGMO.algorithm.mbh.algorithm

      Algorithm to perform mbh 'local' search

---------------

.. class:: PyGMO.algorithm.cstrs_co_evolution

   .. automethod:: PyGMO.algorithm.cstrs_co_evolution.__init__

---------------

.. class:: PyGMO.algorithm.cstrs_immune_system

   .. automethod:: PyGMO.algorithm.cstrs_immune_system.__init__

   .. attribute:: PyGMO.algorithm.cstrs_immune_system.screen_output

      When True, the algorithms produces output on screen 

---------------

.. class:: PyGMO.algorithm.cstrs_core

   .. automethod:: PyGMO.algorithm.cstrs_core.__init__

   .. attribute:: PyGMO.algorithm.cstrs_core.screen_output

      When True, the algorithms produces output on screen 

---------------

.. class:: PyGMO.algorithm.cs

   .. automethod:: PyGMO.algorithm.cs.__init__

---------------

.. class:: PyGMO.algorithm.ihs

   .. automethod:: PyGMO.algorithm.ihs.__init__

---------------

.. class:: PyGMO.algorithm.monte_carlo

   .. automethod:: PyGMO.algorithm.monte_carlo.__init__

---------------

.. class:: PyGMO.algorithm.py_example

   .. automethod:: PyGMO.algorithm.py_example.__init__

---------------

.. class:: PyGMO.algorithm.py_cmaes

   .. automethod:: PyGMO.algorithm.py_cmaes.__init__

---------------

.. class:: PyGMO.algorithm.cmaes

   .. automethod:: PyGMO.algorithm.cmaes.__init__

---------------

.. class:: PyGMO.algorithm.scipy_fmin

   .. automethod:: PyGMO.algorithm.scipy_fmin.__init__

---------------

.. class:: PyGMO.algorithm.scipy_l_bfgs_b

   .. automethod:: PyGMO.algorithm.scipy_l_bfgs_b.__init__

---------------

.. class:: PyGMO.algorithm.scipy_slsqp

   .. automethod:: PyGMO.algorithm.scipy_slsqp.__init__

---------------

.. class:: PyGMO.algorithm.scipy_tnc

   .. automethod:: PyGMO.algorithm.scipy_tnc.__init__

   .. attribute:: PyGMO.algorithm.scipy_tnc.screen_output

      When True, the algorithms produces output on screen 

---------------

.. class:: PyGMO.algorithm.scipy_cobyla

   .. automethod:: PyGMO.algorithm.scipy_cobyla.__init__

   .. attribute:: PyGMO.algorithm.scipy_cobyla.screen_output

      When True, the algorithms produces output on screen 

---------------

.. class:: PyGMO.algorithm.nlopt_cobyla

   .. automethod:: PyGMO.algorithm.nlopt_cobyla.__init__

---------------

.. class:: PyGMO.algorithm.nlopt_bobyqa

   .. automethod:: PyGMO.algorithm.nlopt_bobyqa.__init__

---------------

.. class:: PyGMO.algorithm.nlopt_sbplx

   .. automethod:: PyGMO.algorithm.nlopt_sbplx.__init__

---------------

.. class:: PyGMO.algorithm.nlopt_mma

   .. automethod:: PyGMO.algorithm.nlopt_mma.__init__

---------------

.. class:: PyGMO.algorithm.nlopt_auglag

   .. automethod:: PyGMO.algorithm.nlopt_auglag.__init__

---------------

.. class:: PyGMO.algorithm.nlopt_auglag_eq

   .. automethod:: PyGMO.algorithm.nlopt_auglag_eq.__init__

---------------

.. class:: PyGMO.algorithm.nlopt_slsqp

   .. automethod:: PyGMO.algorithm.nlopt_slsqp.__init__

---------------

.. class:: PyGMO.algorithm.gsl_nm2rand

   .. automethod:: PyGMO.algorithm.gsl_nm2rand.__init__

---------------

.. class:: PyGMO.algorithm.gsl_nm2

   .. automethod:: PyGMO.algorithm.gsl_nm2.__init__

---------------

.. class:: PyGMO.algorithm.gsl_nm

   .. automethod:: PyGMO.algorithm.gsl_nm.__init__

---------------

.. class:: PyGMO.algorithm.gsl_pr

   .. automethod:: PyGMO.algorithm.gsl_pr.__init__

---------------

.. class:: PyGMO.algorithm.gsl_fr

   .. automethod:: PyGMO.algorithm.gsl_fr.__init__

---------------

.. class:: PyGMO.algorithm.gsl_bfgs2

   .. automethod:: PyGMO.algorithm.gsl_bfgs2.__init__

---------------

.. class:: PyGMO.algorithm.gsl_bfgs

   .. automethod:: PyGMO.algorithm.gsl_bfgs.__init__

---------------

.. class:: PyGMO.algorithm.snopt

   .. automethod:: PyGMO.algorithm.snopt.__init__

   .. attribute:: PyGMO.algorithm.snopt.screen_output

      When True, the algorithms produces output on screen 

---------------

.. class:: PyGMO.algorithm.ipopt

   .. automethod:: PyGMO.algorithm.ipopt.__init__

   .. attribute:: PyGMO.algorithm.ipopt.screen_output

      When True, the algorithms produces output on screen 

---------------

.. class:: PyGMO.algorithm.cstrs_self_adaptive

   .. automethod:: PyGMO.algorithm.cstrs_self_adaptive.__init__
