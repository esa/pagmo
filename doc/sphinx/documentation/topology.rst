Topology 
=========

The classes of this submodule are all instances of the same base class used to define the migration paths in the
:class:`PyGMO.archipelago`. When a new :class:`PyGMO.island` is pushed back into 
the :class:`PyGMO.archipelago` the various connections are rewired as to respect the topological properties defined 
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
  archi.topology.draw() #requires networkx to be installed
  

.. class:: PyGMO.topology.ring()

   Ring topology (links go in both directions)

.. class:: PyGMO.topology.one_way_ring()

   Ring topology (links go in one direction)

.. class:: PyGMO.topology.fully_connected()

   Fully connected topology (i.e. all islands are connected bi-directionally)

.. class:: PyGMO.topology.pan()

   Islands are connected in a ring, except island No. 0 that is connected to the ring only receiving migrants. This
   topology was first conceived to have global optimization happening on the outer ring while using some local optimizer
   in the island 0 to refine the solution

.. class:: PyGMO.topology.rim()

   Islands are connected in a ring and they all are unidirectinally connected to the No.0 island. This topology was first
   conceived having in mind to put a local optimizer in the island No. 0

.. class:: PyGMO.topology.hypercube()

   Hypercube topology



