Problems
==========

A Quick Look
------------

Problems in PyGMO are objects, first constructed and then used in conjunction to an algorithm.
The user can implement its own problem directly in Python, in which case he needs to inherit from 
:class:`PyGMO.problem.base` or :class:`PyGMO.problem.base_stochastic` class. You may see 
:ref:`adding_a_new_optimization_problem` or :ref:`adding_a_new_stochastic_optimization_problem` 

Meta-problems
^^^^^^^^^^^^^
================================== ========================================= ===========================================
Common Name                        Name in PyGMO                             Comments
================================== ========================================= ===========================================
Rotated                            :class:`PyGMO.problem.rotated`            
Shifted                            :class:`PyGMO.problem.shifted`            
Normalized                         :class:`PyGMO.problem.normalized`        
Noisy                              :class:`PyGMO.problem.noisy`
Decompose                          :class:`PyGMO.problem.decompose`    
Death-penalty                      :class:`PyGMO.problem.death_penalty`      Minimization assumed.
Constrained to MO                  :class:`PyGMO.problem.con2mo`             
Constrained to Unconstrained       :class:`PyGMO.problem.con2uncon`         
================================== ========================================= ===========================================

Box-Constrained Continuous Single-Objective
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
================================== ========================================= ===========================================
Common Name                        Name in PyGMO                             Comments
================================== ========================================= ===========================================
Ackley                             :class:`PyGMO.problem.ackley`            
Bukin F6                           :class:`PyGMO.problem.bukin`              A difficult bi-dimensional problem
Branin                             :class:`PyGMO.problem.branin`             Bi-dimensional problem
CEC2013                            :class:`PyGMO.problem.cec2013`            28 problems part of CEC2013 Competition
De Jong                            :class:`PyGMO.problem.dejong`
De Jong                            :class:`PyGMO.problem.py_example`         Implemented directly in Python
Griewank                           :class:`PyGMO.problem.griewank`
Himmelblau                         :class:`PyGMO.problem.himmelblau`         Bi-dimensional problem
Lennard-Jones                      :class:`PyGMO.problem.lennard_jones`
Michalewicz                        :class:`PyGMO.problem.michalewicz`
Rosenbrock                         :class:`PyGMO.problem.rosenbrock`
Rastrigin                          :class:`PyGMO.problem.rastrigin`
Schwefel                           :class:`PyGMO.problem.schwefel`
MGA-1DSM (tof encoding)            :class:`PyGMO.problem.mga_1dsm_tof`       Requires the GTOP database option active      
MGA-1DSM (alpha encoding)          :class:`PyGMO.problem.mga_1dsm_alpha`     Requires the GTOP database option active      
Cassini 1                          :class:`PyGMO.problem.cassini_1`          Requires the GTOP database option active
Cassini 2                          :class:`PyGMO.problem.cassini_2`          Requires the GTOP database option active
Rosetta                            :class:`PyGMO.problem.rosetta`            Requires the GTOP database option active
Tandem                             :class:`PyGMO.problem.tandem`             Requires the GTOP database option active
Laplace                            :class:`PyGMO.problem.tandem`             Requires the GTOP database option active
Messenger (Full Problem)           :class:`PyGMO.problem.messenger_full`     Requires the GTOP database option active
GTOC1                              :class:`PyGMO.problem.gtoc_1`             Requires the GTOP database option active
Sagas                              :class:`PyGMO.problem.sagas`              Requires the GTOP database option active
================================== ========================================= ===========================================

Box-Constrained Continuous Multi-Objective
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
================================== ========================================= ===========================================
Common Name                        Name in PyGMO                             Comments
================================== ========================================= ===========================================
Kursawe's study                    :class:`PyGMO.problem.kur`
Fonseca and Fleming's study        :class:`PyGMO.problem.fon`
Poloni's study                     :class:`PyGMO.problem.pol`
Shaffer's study                    :class:`PyGMO.problem.sch`
ZDT                                :class:`PyGMO.problem.zdt`
DTLZ                               :class:`PyGMO.problem.dtlz`
CEC2009 (UF1-UF10)                 :class:`PyGMO.problem.cec2009`            UF problems from CEC2009 Competition.
MGA-1DSM (tof encoding)            :class:`PyGMO.problem.mga_1dsm_tof`       Requires the GTOP database option active      
MGA-1DSM (alpha encoding)          :class:`PyGMO.problem.mga_1dsm_alpha`     Requires the GTOP database option active      
Cassini 1                          :class:`PyGMO.problem.cassini_1`          Requires the GTOP database option active
================================== ========================================= ===========================================

Constrained Continuous Single-Objective
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
================================== ========================================= ===========================================
Common Name                        Name in PyGMO                             Comments
================================== ========================================= ===========================================
CEC2006                            :class:`PyGMO.problem.cec2006`            24 problems part of CEC2006 Competition
Pressure vessel design             :class:`PyGMO.problem.pressure_vessel`    
Welded beam design                 :class:`PyGMO.problem.welded_beam`        
Tension compression string design  :class:`PyGMO.problem.tens_comp_string`    
Luksan Vlcek 1                     :class:`PyGMO.problem.luksan_vlcek_1`
Luksan Vlcek 2                     :class:`PyGMO.problem.luksan_vlcek_2`
Luksan Vlcek 3                     :class:`PyGMO.problem.luksan_vlcek_3`
Planet to Planet LT Transfer       :class:`PyGMO.problem.py_pl2pl`           Requires PyKEP. Implemented in Python
SNOPT Toy-Problem                  :class:`PyGMO.problem.snopt_toyprob`      
GTOC2 (Full Problem)               :class:`PyGMO.problem.gtoc_2`             Requires the GTOP database option active
================================== ========================================= ===========================================

Constrained Continuous Multi-Objective
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
================================== ========================================= ===========================================
Common Name                        Name in PyGMO                             Comments
================================== ========================================= ===========================================
CEC2009 (CF1-CF10)                 :class:`PyGMO.problem.cec2009`            CF problems from CEC2009 Competition.
================================== ========================================= ===========================================

Box-Constrained Integer Single-Objective
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
================================== ========================================= ===========================================
Common Name                        Name in PyGMO                             Comments
================================== ========================================= ===========================================
String Match                       :class:`PyGMO.problem.string_match`
================================== ========================================= ===========================================

Constrained Integer Single-Objective
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
================================== ========================================= ===========================================
Common Name                        Name in PyGMO                             Comments
================================== ========================================= ===========================================
Golomb Ruler                       :class:`PyGMO.problem.golomb_ruler`
Traveling Salesman                 :class:`PyGMO.problem.tsp`
================================== ========================================= ===========================================

Stochastic Objective Function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
================================== =============================================== ===========================================
Common Name                        Name in PyGMO                                   Comments
================================== =============================================== ===========================================
Inventory Problem                  :class:`PyGMO.problem.inventory`
MIT SPHERES                        :class:`PyGMO.problem.mit_spheres`
Noisy De Jong                      :class:`PyGMO.problem.py_example_stochastic`
================================== =============================================== ===========================================



Detailed Documentation
----------------------

.. class:: PyGMO.problem.base

   .. automethod:: PyGMO.problem.base.__init__

   .. autoattribute:: PyGMO.problem.base.best_x

   .. autoattribute:: PyGMO.problem.base.best_f

   .. autoattribute:: PyGMO.problem.base.best_c

   .. autoattribute:: PyGMO.problem.base.dimension

   .. autoattribute:: PyGMO.problem.base.f_dimension

   .. autoattribute:: PyGMO.problem.base.i_dimension

   .. autoattribute:: PyGMO.problem.base.c_dimension

   .. autoattribute:: PyGMO.problem.base.ic_dimension

   .. autoattribute:: PyGMO.problem.base.c_tol
   
   .. method:: _objfun_impl(self, x)
   
      This is a virtual function tham must be re-implemented in the derived class and must
      return a tuple packing as many numbers as the problem objectives (n_obj)
      
   .. method:: _compute_constraints_impl(self, x)
   
      This is a virtual function that can be re-implemented in the derived class (if c_dim>0) and must return a tuple 
      packing as many numbers as the declared dimension of the problem constraints (c_dim). 
      Inequality constarints need to be packed at last.

   .. method:: _compare_fitness_impl(self, f1, f2)
   
      This is a virtual function that can be re-implemented in the derived class and must return a boolean value.
      Return true if f1 Pareto dominate f2, false otherwise. The default implementation will
      assume minimisation for each one of the f components i.e., each pair of corresponding elements
      in f1 and f2 is compared: if all elements in f1 are less or equal to the corresponding element in f2 (and
      at least one is less), true will be returned. Otherwise, false will be returned.     
            
   .. method:: _compare_constraints_impl(self, c1, c2)
   
      This is a virtual function tham can be re-implemented in the derived class (if c_dim>0) and must return a boolean value.
      Return true if c1 is a strictly better constraint vector than c2, false otherwise. 
      Default implementation will return true under the following conditions, tested in order: c1 satisfies more constraints than c2,
      c1 and c2 satisfy the same number of constraints and the L2 norm of the constraint mismatches for c1 is smaller than for c2.
      Otherwise, false will be returned.

   .. method:: _compare_fc_impl(self, f1, c1, f2, c2)
   
      This is a virtual function that can be re-implemented in the derived class (if c_dim>0) and must return a boolean value. 
      By default, the function will perform sanity checks on the input arguments and will
      then call _compare_constraints_impl() if the constraint dimensions is not null, _compare_fitness_impl() otherwise.   
      
   .. automethod:: PyGMO.problem.base.reset_caches

   .. automethod:: PyGMO.problem.base.set_bounds
   
   .. automethod:: PyGMO.problem.base.feasibility_x
   
   .. automethod:: PyGMO.problem.base.feasibility_c

-----------------

   
.. class:: PyGMO.problem.death_penalty

   .. automethod:: PyGMO.problem.death_penalty.__init__

-----------------
   
.. class:: PyGMO.problem.con2mo

   .. automethod:: PyGMO.problem.con2mo.__init__

-----------------

.. class:: PyGMO.problem.con2uncon

   .. automethod:: PyGMO.problem.con2uncon.__init__

-----------------

.. class:: PyGMO.problem.shifted

   .. automethod:: PyGMO.problem.shifted.__init__
   
   .. attribute:: shift_vector
   
      The shift vector defining the new problem
      
   .. method:: PyGMO.problem.shifted.deshift((tuple) x)

      Returns the de-shifted decision vector

-----------------
   
.. class:: PyGMO.problem.rotated

   .. automethod:: PyGMO.problem.rotated.__init__
   
   .. attribute:: rotation
   
      The rotation matrix defining the new problem
      
   .. method:: PyGMO.problem.rotated.derotate((tuple) x)

      Returns the de-rotated decision vector

-----------------

      
.. class:: PyGMO.problem.noisy

   .. automethod:: PyGMO.problem.noisy.__init__

-----------------

.. class:: PyGMO.problem.normalized

   .. automethod:: PyGMO.problem.normalized.__init__

   .. method:: PyGMO.problem.normalized.denormalize((tuple) x)

      Returns the de-normalized decision vector

-----------------

.. class:: PyGMO.problem.decompose

   .. automethod:: PyGMO.problem.decompose.__init__
   
   .. attribute:: weights
      
      The weights vector

   .. attribute:: PyGMO.problem.decompose.WEIGHTED
      
      Weighted decomposition method

   .. attribute:: PyGMO.problem.decompose.TCHEBYCHEFF
      
      Tchebycheff decomposition method

   .. attribute:: PyGMO.problem.decompose.BI
      
      Boundary Intersection decomposition method

-----------------
   
.. class:: PyGMO.problem.ackley

   .. automethod:: PyGMO.problem.ackley.__init__

-----------------
   
.. class:: PyGMO.problem.bukin

   .. automethod:: PyGMO.problem.bukin.__init__

-----------------

.. class:: PyGMO.problem.cec2006

   .. automethod:: PyGMO.problem.cec2006.__init__

-----------------

.. class:: PyGMO.problem.pressure_vessel

   .. automethod:: PyGMO.problem.pressure_vessel.__init__

-----------------

.. class:: PyGMO.problem.welded_beam

   .. automethod:: PyGMO.problem.welded_beam.__init__

-----------------

.. class:: PyGMO.problem.tens_comp_string

   .. automethod:: PyGMO.problem.tens_comp_string.__init__

-----------------

.. class:: PyGMO.problem.cec2009

   .. automethod:: PyGMO.problem.cec2009.__init__

-----------------
   
.. class:: PyGMO.problem.cec2013

   .. automethod:: PyGMO.problem.cec2013.__init__

-----------------

.. class:: PyGMO.problem.rosenbrock

   .. automethod:: PyGMO.problem.rosenbrock.__init__

-----------------

.. class:: PyGMO.problem.string_match

   .. automethod:: PyGMO.problem.string_match.__init__

   .. method:: PyGMO.problem.pretty(x)
 
      Returns a string decoding the chromosome

-----------------

.. class:: PyGMO.problem.rastrigin

   .. automethod:: PyGMO.problem.rastrigin.__init__

-----------------

.. class:: PyGMO.problem.schwefel

   .. automethod:: PyGMO.problem.schwefel.__init__

-----------------

.. class:: PyGMO.problem.dejong

   .. automethod:: PyGMO.problem.dejong.__init__

-----------------

.. class:: PyGMO.problem.py_example

   .. automethod:: PyGMO.problem.py_example.__init__

-----------------

.. class:: PyGMO.problem.griewank

   .. automethod:: PyGMO.problem.griewank.__init__

-----------------

.. class:: PyGMO.problem.lennard_jones

   .. automethod:: PyGMO.problem.lennard_jones.__init__

-----------------

.. class:: PyGMO.problem.branin

   .. automethod:: PyGMO.problem.branin.__init__

-----------------

.. class:: PyGMO.problem.himmelblau

   .. automethod:: PyGMO.problem.himmelblau.__init__

-----------------

.. class:: PyGMO.problem.michalewicz

   .. automethod:: PyGMO.problem.michalewicz.__init__

-----------------

.. class:: PyGMO.problem.kur

   .. automethod:: PyGMO.problem.kur.__init__

-----------------

.. class:: PyGMO.problem.fon

   .. automethod:: PyGMO.problem.fon.__init__

-----------------

.. class:: PyGMO.problem.pol

   .. automethod:: PyGMO.problem.pol.__init__

-----------------

.. class:: PyGMO.problem.sch

   .. automethod:: PyGMO.problem.sch.__init__

-----------------

.. class:: PyGMO.problem.zdt

   .. automethod:: PyGMO.problem.zdt.__init__

   .. automethod:: PyGMO.problem.zdt.p_distance

-----------------
  
.. class:: PyGMO.problem.dtlz

   .. automethod:: PyGMO.problem.dtlz.__init__
   
   .. automethod:: PyGMO.problem.dtlz.p_distance
   
   .. automethod:: PyGMO.problem.dtlz.plot

-----------------

.. class:: PyGMO.problem.tsp

   .. automethod:: PyGMO.problem.tsp.__init__

-----------------

.. class:: PyGMO.problem.golomb_ruler

   .. automethod:: PyGMO.problem.golomb_ruler.__init__

-----------------

.. class:: PyGMO.problem.luksan_vlcek_1

   .. automethod:: PyGMO.problem.luksan_vlcek_1.__init__

-----------------

.. class:: PyGMO.problem.luksan_vlcek_2

   .. automethod:: PyGMO.problem.luksan_vlcek_2.__init__

-----------------

.. class:: PyGMO.problem.luksan_vlcek_3

   .. automethod:: PyGMO.problem.luksan_vlcek_3.__init__

-----------------

.. class:: PyGMO.problem.snopt_toyprob

   .. automethod:: PyGMO.problem.snopt_toyprob.__init__

-----------------

.. class:: PyGMO.problem.inventory

   .. automethod:: PyGMO.problem.inventory.__init__

-----------------

.. class:: PyGMO.problem.py_example_stochastic

   .. automethod:: PyGMO.problem.py_example_stochastic.__init__

-----------------

.. class:: PyGMO.problem.mit_spheres

   .. automethod:: PyGMO.problem.mit_spheres.__init__

   .. method:: post_evaluate((tuple) x, (int) N, (int) seed) -> (tuple) out

      Returns a tuple with the N post evaluation results of chromosome x w.r.t. conditions generated by seed.
      The returned tuple has the structure [ic, fit] and is sorted by fit. Where ic are the initial conditions and fit the
      Evaluated fitness. 

   .. method:: simulate((tuple) x, (tuple) ic, (int) N) -> (tuple) world_states

      Returns the SPHERES coordinates as evaluated in one simulation with initial conditions ic and in
      N points 

   .. method:: visualize((tuple) world_states)

      Requires VPython installed. It opens a graphical display and animate the motion of the three SPHERES
      as desribed by the world_state tuple (output from the simulate method)

-----------------

.. class:: PyGMO.problem.mga_1dsm_tof

   .. method:: PyGMO.problem.mga_1dsm_tof.__init__(seq, t0, tof, vinf, multi_objective=False, add_vinf_dep=False, add_vinf_arr=True)
   
    Constructs an mga_1dsm problem (tof-encoding)

    * seq: list of PyKEP planets defining the encounter sequence, including the starting planet (default: earth venus earth)
    * t0: list of two epochs defining the launch window (default: 2000-Jan-01 00:00:00 to 2002-Sep-27 00:00:00)
    * tof: list of intervals defining the times of flight in days (default: [[50,900],[50,900]])
    * vinf: list of two floats defining the minimum and maximum allowed initial hyperbolic velocity at launch in km/sec (default: [0.5, 2.5])
    * multi_objective: when True constructs a multiobjective problem (dv, T)
    * add_vinf_dep: when True the computed Dv includes the initial hyperbolic velocity (at launch)
    * add_vinf_arr: when True the computed Dv includes the final hyperbolic velocity (at arrival)

    USAGE: problem.mga_1dsm(seq = [planet_ss('earth'),planet_ss('venus'),planet_ss('earth')], t0 = [epoch(0),epoch(1000)], tof = [ [200, 700], [200, 700] ], vinf = [0.5, 2.5], multi_objective = False, add_vinf_dep = False, add_vinf_arr = True)

    .. automethod:: PyGMO.problem.mga_1dsm_tof.plot

    .. method:: PyGMO.problem.mga_1dsm_tof.set_tof(tof)

    Resets the tof-bounds by the provided list of epochs. Needs a list consisting of lower/upper bound tuples.

    .. method:: PyGMO.problem.mga_1dsm_tof.set_launch_window((tuple) t0)

    Resets the launch windows to the lower and upper bounds given by tuple t0. Bounds need to be epochs.

    .. method:: PyGMO.problem.mga_1dsm_tof.set_vinf((double) vinf_u)

    Sets the upper bound for vinf to vinf_u

    .. method:: PyGMO.problem.mga_1dsm_tof.pretty((tuple) x) -> (string) out

    Returns a string with informations about tour encoded by x

-----------------

.. class:: PyGMO.problem.mga_1dsm_alpha

   .. method:: PyGMO.problem.mga_1dsm_alpha.__init__(seq, t0, tof, vinf, multi_objective=False, add_vinf_dep=False, add_vinf_arr=True)
   
    Constructs an mga_1dsm problem (alpha-encoding)

    * seq: list of PyKEP planets defining the encounter sequence, including the starting planet (default: earth venus earth)
    * t0: list of two epochs defining the launch window (default: 2000-Jan-01 00:00:00 to 2002-Sep-27 00:00:00)
    * tof: list of two floats defining the minimum and maximum allowed mission length in days (default: [365.25, 1826.35])
    * vinf: list of two floats defining the minimum and maximum allowed initial hyperbolic velocity at launch in km/sec (default: [0.5, 2.5])
    * multi_objective: when True constructs a multiobjective problem (dv, T)
    * add_vinf_dep: when True the computed Dv includes the initial hyperbolic velocity (at launch)
    * add_vinf_arr: when True the computed Dv includes the final hyperbolic velocity (at arrival)

    USAGE: problem.mga_1dsm(seq = [planet_ss('earth'),planet_ss('venus'),planet_ss('earth')], t0 = [epoch(0),epoch(1000)], tof = [ [200, 700], [200, 700] ], vinf = [0.5, 2.5], multi_objective = False, add_vinf_dep = False, add_vinf_arr = True)

    .. automethod:: PyGMO.problem.mga_1dsm_alpha.plot

    .. method:: PyGMO.problem.mga_1dsm_alpha.set_tof((tuple) tof)

    Resets the tof-bounds by the provided tuple of epochs.

    .. method:: PyGMO.problem.mga_1dsm_alpha.set_launch_window((tuple) t0)

    Resets the launch windows to the lower and upper bounds given by tuple t0. Bounds need to be epochs.

    .. method:: PyGMO.problem.mga_1dsm_alpha.set_vinf((double) vinf_u)

    Sets the upper bound for vinf to vinf_u

    .. method:: PyGMO.problem.mga_1dsm_alpha.pretty((tuple) x) -> (string) out

    Returns a string with informations about tour encoded by x

-----------------

.. class:: PyGMO.problem.cassini_1

   .. automethod:: PyGMO.problem.cassini_1.__init__

-----------------

.. class:: PyGMO.problem.cassini_2

   .. automethod:: PyGMO.problem.cassini_2.__init__

-----------------

.. class:: PyGMO.problem.messenger_full

   .. automethod:: PyGMO.problem.messenger_full.__init__

-----------------

.. class:: PyGMO.problem.rosetta

   .. automethod:: PyGMO.problem.rosetta.__init__

-----------------

.. class:: PyGMO.problem.laplace

   .. automethod:: PyGMO.problem.laplace.__init__

-----------------

.. class:: PyGMO.problem.tandem

   .. automethod:: PyGMO.problem.tandem.__init__

-----------------

.. class:: PyGMO.problem.gtoc_1

   .. automethod:: PyGMO.problem.gtoc_1.__init__

-----------------

.. class:: PyGMO.problem.gtoc_2

   .. automethod:: PyGMO.problem.gtoc_2.__init__

-----------------

.. class:: PyGMO.problem.py_pl2pl

   .. automethod:: PyGMO.problem.py_pl2pl.__init__

-----------------

.. class:: PyGMO.problem.sagas

   .. automethod:: PyGMO.problem.sagas.__init__
