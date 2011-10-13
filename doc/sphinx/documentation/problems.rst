Problems
==========

A Quick Look
------------

Problems in PyGMO are objects, first constructed and then used in conjunction to an algorithm. The user can implement its own problem directly
in Python, in which case he needs to inherit from the :class:`PyGMO.problem.base` class. 


Box-Constrained Continuous Single-Objective
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
================================== ========================================= ===========================================
Common Name                        Name in PyGMO                             Comments
================================== ========================================= ===========================================
Ackley                             :class:`PyGMO.problem.ackley`                 
Rosenbrock                         :class:`PyGMO.problem.rosenbrock`
Rastrigin                          :class:`PyGMO.problem.rastrigin`
Schwefel                           :class:`PyGMO.problem.schwefel`
De Jong                            :class:`PyGMO.problem.dejong`
De Jong                            :class:`PyGMO.problem.py_example`         Implemented directly in Python
Griewank                           :class:`PyGMO.problem.griewank`
Lennard-Jones                      :class:`PyGMO.problem.lennard_jones`
Branin                             :class:`PyGMO.problem.branin`             Bi-dimensional problem
Himmelblau                         :class:`PyGMO.problem.himmelblau`         Bi-dimensional problem
Michalewicz                        :class:`PyGMO.problem.michalewicz`
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
ZDT1                               :class:`PyGMO.problem.zdt1`
ZDT2                               :class:`PyGMO.problem.zdt2`         
ZDT4                               :class:`PyGMO.problem.zdt4`
ZDT6                               :class:`PyGMO.problem.zdt6`
================================== ========================================= ===========================================

Constrained Continuous Single-Objective
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
================================== ========================================= ===========================================
Common Name                        Name in PyGMO                             Comments
================================== ========================================= ===========================================
Luksan Vlcek 1                     :class:`PyGMO.problem.luksan_vlcek_1`
Luksan Vlcek 2                     :class:`PyGMO.problem.luksan_vlcek_2`
Luksan Vlcek 3                     :class:`PyGMO.problem.luksan_vlcek_3`
Planet to Planet LT Transfer       :class:`PyGMO.problem.py_pl2pl`           Requires PyKEP. Implemented in Python
SNOPT Toy-Problem                  :class:`PyGMO.problem.snopt_toyprob`      
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
Knapsack                           :class:`PyGMO.problem.knapsack`
================================== ========================================= ===========================================

Stochastic Objective Function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
================================== ========================================= ===========================================
Common Name                        Name in PyGMO                             Comments
================================== ========================================= ===========================================
Inventory Problem                  :class:`PyGMO.problem.inventory`
mit_spheres                        :class:`PyGMO.problem.mit_spheres`
================================== ========================================= ===========================================



Detailed Documentation
----------------------

.. autoclass:: PyGMO.problem.ackley

   .. automethod:: PyGMO.problem.ackley.__init__

.. autoclass:: PyGMO.problem.rosenbrock

   .. automethod:: PyGMO.problem.rosenbrock.__init__

.. autoclass:: PyGMO.problem.string_match

   .. automethod:: PyGMO.problem.string_match.__init__

   .. method:: PyGMO.problem.pretty(x)
 
      Returns a string decoding the chromosome

.. autoclass:: PyGMO.problem.rastrigin

   .. automethod:: PyGMO.problem.rastrigin.__init__

.. autoclass:: PyGMO.problem.schwefel

   .. automethod:: PyGMO.problem.schwefel.__init__

.. autoclass:: PyGMO.problem.dejong

   .. automethod:: PyGMO.problem.dejong.__init__

.. autoclass:: PyGMO.problem.py_example

   .. automethod:: PyGMO.problem.py_example.__init__

.. autoclass:: PyGMO.problem.griewank

   .. automethod:: PyGMO.problem.griewank.__init__

.. autoclass:: PyGMO.problem.lennard_jones

   .. automethod:: PyGMO.problem.lennard_jones.__init__

.. autoclass:: PyGMO.problem.branin

   .. automethod:: PyGMO.problem.branin.__init__

.. autoclass:: PyGMO.problem.himmelblau

   .. automethod:: PyGMO.problem.himmelblau.__init__

.. autoclass:: PyGMO.problem.michalewicz

   .. automethod:: PyGMO.problem.michalewicz.__init__

.. autoclass:: PyGMO.problem.kur

   .. automethod:: PyGMO.problem.kur.__init__

.. autoclass:: PyGMO.problem.fon

   .. automethod:: PyGMO.problem.fon.__init__

.. autoclass:: PyGMO.problem.pol

   .. automethod:: PyGMO.problem.pol.__init__

.. autoclass:: PyGMO.problem.sch

   .. automethod:: PyGMO.problem.sch.__init__

.. autoclass:: PyGMO.problem.zdt1

   .. automethod:: PyGMO.problem.zdt1.__init__

.. autoclass:: PyGMO.problem.zdt2

   .. automethod:: PyGMO.problem.zdt2.__init__

.. autoclass:: PyGMO.problem.zdt4

   .. automethod:: PyGMO.problem.zdt4.__init__

.. autoclass:: PyGMO.problem.zdt6

   .. automethod:: PyGMO.problem.zdt6.__init__

.. autoclass:: PyGMO.problem.tsp

   .. automethod:: PyGMO.problem.tsp.__init__

.. autoclass:: PyGMO.problem.golomb_ruler

   .. automethod:: PyGMO.problem.golomb_ruler.__init__

.. autoclass:: PyGMO.problem.knapsack

   .. automethod:: PyGMO.problem.knapsack.__init__

.. autoclass:: PyGMO.problem.luksan_vlcek_1

   .. automethod:: PyGMO.problem.luksan_vlcek_1.__init__

.. autoclass:: PyGMO.problem.luksan_vlcek_2

   .. automethod:: PyGMO.problem.luksan_vlcek_2.__init__

.. autoclass:: PyGMO.problem.luksan_vlcek_3

   .. automethod:: PyGMO.problem.luksan_vlcek_3.__init__

.. autoclass:: PyGMO.problem.snopt_toyprob

   .. automethod:: PyGMO.problem.snopt_toyprob.__init__

.. autoclass:: PyGMO.problem.inventory

   .. automethod:: PyGMO.problem.inventory.__init__

.. autoclass:: PyGMO.problem.mit_spheres

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

      