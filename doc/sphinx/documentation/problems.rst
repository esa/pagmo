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
