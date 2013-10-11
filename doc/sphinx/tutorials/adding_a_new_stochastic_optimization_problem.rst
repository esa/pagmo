.. _adding_a_new_stochastic_optimization_problem:

================================================================
Adding a new stochastic optimization problem
================================================================

Assuming you have made the :ref:`adding_a_new_optimization_problem` tutorial, 
this second tutorial should be rather straight forward. We will learn how to code
a stochastic optimization problem, that is a problem where the objective function is stochastic.
In other words, the objective function depends on pseudo-random numbers.

In a nutshell .... we will write a class deriving from problem.base_stochastic
and reimplement some of its 'virtual' methods. With respect to problem.base, this base
class has an added data member called self.seed which can be used to control the pseudo random numbers
and that is changed by the algorithms that are compatible with stochastic optimization problems.

.. code-block:: python

   from PyGMO.problem import base_stochastic
   class my_problem_stochastic(base_stochastic):
	"""
	Noisy De Jong (sphere) function implemented purely in Python.
	
	USAGE: my_problem_stochastic(dim = 10, seed=0)

	* dim problem dimension
	* seed initial random seed
	"""
	def __init__(self, dim = 10, seed = 0):
		#First we call the constructor of the base stochastic class. (Only
		#unconstrained single objective problems can be stochastic in PyGMO)
		super(my_problem_stochastic,self).__init__(dim, seed)

		#then we set the problem bounds (in this case equal for all components)
		self.set_bounds(-5.12,5.12)

		#and we define some additional 'private' data members (not really necessary in
		#this case, but ... hey this is a tutorial)
		self.__dim = dim

	def _objfun_impl(self,x):
		from random import random as drng
		from random import seed

		#We initialize the random number generator using the 
		#data member seed (in base_stochastic). This will be changed by suitable
		#algorithms when a stochastic problem is used. The mod operation avoids overflows
		
		seed(self.seed)
		
		#We write the objfun using the same pseudorandonm sequence
		#as long as self.seed is unchanged.
		f = 0;
		for i in range(self.__dim):
			noise = (2 * drng() - 1) / 10
			f = f + (x[i] + noise)*(x[i] + noise)
		return (f,)
	def human_readable_extra(self):
		return "\n\tSeed: " + str(self.seed)

We may then put this in a file, say my_module.py and use Particle Swarm Optimization .... 
the generational version ... to solve it

.. code-block:: python

   from PyGMO import *
   from  my_module import my_problem_stochastic

   prob = my_problem_stochastic(dim=10, seed=123456)
   algo = algorithm.pso_gen(gen=1)
   isl = island(algo,prob,20)
   for i in range(30):
	isl.evolve(1)
	isl.population.champion.f

.. code-block:: python

    (7.636346361215645,)
    (2.9207347362223715,)
    (1.0075035416057239,)
    (0.3461345536433724,)
    (0.148618031035022,)
    (0.08653472558404088,)
    (0.048504492499211634,)
    (0.017069267823961475,)
    (0.032427061740872294,)
    (0.018818646172131907,)
    (0.025077134049593254,)

You can also check that the problem seed has actually changed (the algo does this) by printing 
the seed to screen

.. code-block:: python

   print prob.seed
