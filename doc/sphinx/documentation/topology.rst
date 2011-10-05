Topology 
=========

The classes of this submodule are all instances of the same base class used to define the migration paths in the
:class:`PyGMO.archipelago`. When a new :class:`PyGMO.island` is pushed back into 
the :class:`PyGMO.archipelago` the various connections are rewired as to respect the topological properties defined 
by these  classes. In the picture below a :class:`PyGMO.archipelago` is shown with all its :class:`PyGMO.island` and
the migration paths among them. 

The user will mainly use these classes in the following way:

.. code-block:: python

  from PyGMO import *
  prob = problem.lennard_jones(5)
  algo = algorithm.de(gen = 10) #instantiates artificial bee colony with default params and 10 generations
  topo = topology.ring() #defines a ring topology
  archi = archipelago(algo,prob,8,20,topology=topo) #connects 8 islands in a ring topology
  archi.topology.draw() #requires networkx to be installed
  
.. class:: PyGMO.topology.unconnected()

   Unconnected topology (this corresponds to do parallel independent runs of the algorithms)
  
   .. image:: ../images/unconnected.png

   .. method:: PyGMO.topology.unconnected.draw()

      Uses Python networx module to draw the topology in a graphic window

.. class:: PyGMO.topology.ring()

   Ring topology (links go in both directions)
  
   .. image:: ../images/ring.png

   .. method:: PyGMO.topology.ring.draw()

      Uses Python networx module to draw the topology in a graphic window

.. class:: PyGMO.topology.one_way_ring()

   Ring topology (links go only in one direction)

   .. image:: ../images/one_way_ring.png

   .. method:: PyGMO.topology.one_way_ring.draw()

      Uses Python networx module to draw the topology in a graphic window

.. class:: PyGMO.topology.fully_connected()

   Fully connected topology (i.e. all islands are connected bi-directionally)

   .. image:: ../images/fully_connected.png

   .. method:: PyGMO.topology.fully_connected.draw()

      Uses Python networx module to draw the topology in a graphic window

.. class:: PyGMO.topology.pan

   Islands are connected in a ring, except island No. 0 that is connected to the ring only receiving migrants. This
   topology was first conceived to have global optimization happening on the outer ring while using some local optimizer
   in the island 0 to refine the solution

   .. image:: ../images/pan.png

   .. method:: PyGMO.topology.pan.draw()

      Uses Python networx module to draw the topology in a graphic window

.. class:: PyGMO.topology.rim

   Islands are connected in a ring and they all are also connected to the No.0 island.

   .. image:: ../images/rim.png

   .. method:: PyGMO.topology.rim.draw()

      Uses Python networx module to draw the topology in a graphic window

.. class:: PyGMO.topology.hypercube()

   Hypercube topology

   .. image:: ../images/hypercube.png

   .. method:: PyGMO.topology.hypercube.draw()

      Uses Python networx module to draw the topology in a graphic window

.. class:: PyGMO.topology.watts_strogatz

   The Watts-Strogatz topology is a ring lattice network in which forward
   edges are rewired with random probability. Such a network has small-world properties, 
   including short average path
   lengths and high clustering. When the push_back method is used all the connections 
   are rewired. 

   `Watts-Strogatz topology in wikipedia <http://en.wikipedia.org/wiki/Watts_and_Strogatz_model>`_

   .. method:: PyGMO.topology.watts_strogatz.__init__([(int) K=10, (double)beta=0.1, (int)size=0])

      Builds a Watts_strogatz topology with K neighbours (K/2 on each side) in which forward
      edges are rewired with random probability beta. Since the addition of a single element to the topology implies the
      rewiring of the whole topology, for archipelago objects of large size it is advisable
      to build a topology object outside the archipelago and then assign it to the archipelago

   .. image:: ../images/watts_strogatz.png

   .. method:: PyGMO.topology.watts_strogatz.draw()

      Uses Python networx module to draw the topology in a graphic window

.. class:: PyGMO.topology.erdos_reny

   A random graph

   .. method:: PyGMO.topology.erdos_reny.__init__([(double)p = 0.01)

      Builds a random graph with a probability of p for each bidirectional link

   .. image:: ../images/erdos.png

   .. method:: PyGMO.topology.erdos_reny.draw()

      Uses Python networx module to draw the topology in a graphic window

.. class:: PyGMO.topology.barabasi_albert([(int)m0=3, (int)m=2])

   Topology based on the Barabási-Albert (BA) model

   .. method:: PyGMO.topology.barabasi_albert.__init__([(int)m0=3, (int)m=2])

      Constructs a Barabasi-Albert topology. The construction consists internally of two phases:

      * The first m0 elements added to the network constitute a kernel of nodes connected to each other
        with high probability;
 
      * After the kernel is built, the next elements added to the network
        are connected randomly to m of the existing nodes; the probability of connection 
        is biased linearly towards the most connected nodes.

   .. image:: ../images/ba.png

   .. method:: PyGMO.topology.barabasi_albert.draw()

      Uses Python networx module to draw the topology in a graphic window

.. class:: PyGMO.topology.custom()

   A custom topology. Allows for the construction of any topology via its unique methods.

   .. method:: PyGMO.topology.custom.push_back()
    
      Adds a vertex

   .. method:: PyGMO.topology.custom.add_edge((int)id1, (int)id2)

      Adds a directed adge from vertex id1 to vertex id2

   .. method:: PyGMO.topology.custom.remove_edge((int)id1, (int)id2)

      Removes the directed adge from vertex id1 to vertex id2

   .. method:: PyGMO.topology.custom.remove_all_edges()

      Removes all_edges

   .. method:: PyGMO.topology.custom.draw()

      Uses Python networx module to draw the topology in a graphic window

   Example: 

   .. code-block:: python

      from PyGMO import *
      prob = problem.lennard_jones(5)
      algo = algorithm.de(gen = 10)     #instantiates artificial bee colony with default params and 10 generations
      topo = topology.custom()    #builds a custom topology
      for i in range(30):
           topo.push_back()       #Now we have an unconnected topology of 30 vertices
      topo.add_edge(1,2)
      ....
      topo.add_edge(22,0)
      archi = archipelago(algo,prob,30,20) #constructs an archipelago (we cannot assign here directly topo 
                                           #to the kwarg topology as the archipelago constructor only takes empty topologies
      archi.topology = topo                #sets the topology to the customly constructed one
      archi.topology.draw()                #Draws the topology (this requires networkx to be installed)

