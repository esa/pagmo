Topology 
=========

The classes of this submodule are all instances of the same base class used to define the migration paths in the
:class:`PyGMO.archipelago`. When a new :class:`PyGMO.island` is pushed back into 
the :class:PyGMO.archipelago` the various connections are rewired as to respect the topological properties defined 
by these  classes. In the picture below a :class:`PyGMO.archipelago` is shown with all its :class:`PyGMO.island` and
the migration paths among them. 

.. image:: ../images/BA_topology.png

The user will mainly use these classes in the following way:

.. code-block:: python

  from PyGMO import *
  prob = problem.lennard_jones(5)
  algo = algorithm.de(10) #instantiates artificial bee colony with default params and 10 generations
  topo = topology.ring() #defines a ring topology
  archi = archipelago(algo,prob,8,20,topology=topo) #connects 8 islands in a ring topology
  archi.topology.draw()
  

.. class:: PyGMO.topology.ring()

   Ring topology (links go in both directions)

.. class:: PyGMO.topology.one_way_ring()

   Ring topology (links go in one directions)
