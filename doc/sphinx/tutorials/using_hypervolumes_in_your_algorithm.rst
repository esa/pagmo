.. _using_hypervolumes_in_your_algorithm:

================================================================
Using hypervolumes in your own algorithm
================================================================

In this tutorial we will show how to use the hypervolume features in your own custom evolutionary algorithms.
If you haven't already, we suggest getting familiar with the tutorial :ref:`adding_a_new_algorithm`.

As it is often the case in the tutorials, we will present a fairly simple MOO algorithm, which we would not expect to be particularly good.

.. code-block:: python

  class my_hv_moo_alg(algorithm.base):
  	"""
  	A custom steady-state algorithm, based on the hypervolume computation.
  	"""
  
  	def __init__(self, gen = 10):
  	   """
  	   Constructs a Monte-Carlo (random sampling) algorithm
  
  	   USAGE: my_hv_moo_alg(gen = 10)
  
  	   NOTE: 1. Generates a new individual. 2.Selects the least contributor of the population, and removes it.
  
  	   * gen: number of generations
  	   """
  	   #We start calling the base constructor
  	   super(my_hv_moo_alg,self).__init__()
  	   # Store the number of generations
  	   self.__gen = gen
  
  	# Performs a very simple crossover step
  	def cross(self, ind1, ind2):
  		import random
  		x1 = ind1.cur_x
  		x2 = ind2.cur_x
  		return tuple(random.choice((x1[i], x2[i],)) for i in xrange(len(x1)))
  
  
  	# Evolve method
  	def evolve(self,pop):
  		#If the population is empty (i.e. no individuals) nothing happens
  		if len(pop) == 0:
  		 	 return pop
  
  		import random
  		#The algorithm now starts manipulating the population
  		for s in range(self.__gen):
  			# Initiate new individual by a crossover of two random individuals
  			idx1 = random.randint(0, len(pop) - 1)
  			idx2 = (idx1 + random.randint(0, len(pop) - 2)) % len(pop)
  			ind1 = pop[idx1]
  			ind2 = pop[idx2]
  
  			new_x = self.cross(ind1, ind2)
  			pop.push_back(new_x)
  
  			# Remove the least contributor
  			hv = hypervolume(pop)
  			ref_point = hv.get_nadir_point(1.0)
  			lc_idx = hv.least_contributor(ref_point)
  
  			pop.erase(lc_idx)
  		return pop
  
  	def get_name(self):
  		   return "Custom HV-based MOO"

  prob = problem.dtlz3(fdim=3)
  alg = my_hv_moo_alg(gen = 100)
  pop = population(prob, 100)
  for _ in xrange(100):
    pop = alg.evolve(pop)
    print prob.p_distance(pop)

The algorithm does the following in the *evolve* method:

#. Establish a new individual by performing a very simple crossover on two random individuals
#. Push it to the population
#. Establish the least contributor using the `PyGMO.hypervolume` module
#. Remove the least contributor from the population

Script above should produce an output similar to the one below:

.. code-block:: bash

  P-Distance: 982.90490, Hypervolume: 26703057342.25552
  P-Distance: 923.07177, Hypervolume: 26765472304.25195
  .
  .
  .
  P-Distance: 80.06343, Hypervolume: 26998134491.62578
  P-Distance: 79.00453, Hypervolume: 26998134491.62579

The end effect is not spectacular as the algorithm itself is terribly simple, yet we can observe an improvement over the consecutive generations both in the P-Distance and the hypervolume.
