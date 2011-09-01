Individual and Champion
=======================

.. class:: PyGMO.individual()

   This class represents the concept of an *individual*. In PyGMO, an *individual* is the solution to some optimization
   problem, enriched with some memory of its past (this usually records the best position so far occupied in the solution
   space by a certain individual), and a velocity, representing the variation in the search space of an *individual* position.

   While it is important to know what an *individual* is, the user rarely constructs or manipulate individuals as these
   actions are performed by other classes that also ensure the information contained in an *individual* is consistent
   (i.e. length of the chromosomes, etc. )

   .. method:: __init__()

      Constructs an empty *individual* 

   .. attribute:: cur_f

      Tuple containing the current fitness of the *individual* w.r.t a problem (can contain more than one fitness score for multi-objective optimization)

   .. attribute:: cur_x
      
      Tuple containing the current chromosome (or decision vector) defining the *individual*

   .. attribute:: cur_c

      Tuple containing the current constraint vector defining an *individual* w.r.t a problem (this will be empty in unconstrained optimizaion problems)

   .. attribute:: best_f

      Tuple containing the *individual* fitness corresponding to best_x (the individual can move out of good points as a result of the optimization process)

   .. attribute:: best_x

      Tuple storing the chromosome corresponding to the best position encountered so far since the individual creation (individual can move out of good areas as a resualt of the optimization process)

   .. attribute:: best_c

      Tuple containing the constraint vector corresponding to best_x 

   .. attribute:: cur_v

      Velocity of an *individual* (that is the difference in cur_x between generations). This information is crucial in algorithms such as Particle Swarm Optimization (PSO), but it gets, by default, updated also by other algorithms.


.. class:: PyGMO.champion()

   This class represents the concept of a *champion*, i.e. the best among a set of individuals. Differently from
   an :class:`PyGMO.individual` a *champion* does not have memory, nor a velocity. Similarly from the class :class:`PyGMO.individual` the user rarely constructs or manipulate 
   objects from this class as these actions are performed by other classes that also ensure the information 
   contained in a *champion* is consistent (i.e. it actually is the best in a set of individuals, etc. )

   .. method:: __init__()

      Constructs an empty *champion* 

   .. attribute:: f

      Tuple containing the fitness of the *champion* w.r.t a problem (can contain more than one fitness score for multi-objective optimization)

   .. attribute:: x
      
      Tuple containing the chromosome (or decision vector) defining the *champion*

   .. attribute:: c

      Tuple containing the constraint vector defining a *champion* w.r.t a problem (this will be empty in unconstrained optimizaion problems)
  



