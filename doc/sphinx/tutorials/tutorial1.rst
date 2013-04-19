.. _tutorial1:

================================================================
Tutorial 1: Coding an optimization problem in Python
================================================================

In this Tutorial we will learn how to code a simple optimization problem
(continuous, single objective, unconstrained), so that PyGMO can then apply all of its
algorithmic power to solve it.

In a nutshell .... we will write a class deriving from problem.base
and reimplement some of its 'virtual' methods.

.. code-block:: python

   from PyGMO.problem import base
   class my_problem(base):
	"""
	De Jong (sphere) function implemented purely in Python.
	
	USAGE: my_problem(dim = 10)

	* dim problem dimension
	"""
	def __init__(self, dim = 10):
		#First we call the constructor of the base class telling
		#essentially to PyGMO what kind of problem to expect (1 objective, 0 contraints etc.)
		super(my_problem,self).__init__(dim)

		#then we set the problem bounds (in this case equal for all components)
		self.set_bounds(-5.12,5.12)
		
		#and we define some additional 'private' data members (not really necessary in
		#this case, but ... hey this is a tutorial)
		self.__dim = dim

	#We reimplement the virtual method that defines the objective function.
	def _objfun_impl(self,x):
		f = 0;
		for i in range(self.__dim):
			f = f + (x[i])*(x[i])
		#note that we return a tuple with one element only. In PyGMO the objective functions
		#return tuples so that multi-objective optimization is also possible.
		return (f,)

	#Finally we also reimplement a virtual method that adds some output to the __repr__ method
	def human_readable_extra(self):
		return "\n\t Problem dimension: " + str(self.__dim)

Note that by default PyGMO will assume one wants to minimize the objective function. 

We may then put this in a file, say my_module.py and use, for example, Artificial Bee Colony .... with
20 individuals ....

.. code-block:: python

   from PyGMO import *
   from  my_module import my_problem

   prob = my_problem(dim=10)
   algo = algorithm.bee_colony(gen=100)
   isl = island(algo,prob,20)
   isl.evolve(1); isl.join()
   print isl.population.champion.f

And we are done!!!!! (the output will be something like 10^-27, no big deal for a sphere problem)

NOTE1: This simple tutorial is implemented in PyGMO under the name PyGMO.problem.py_example
NOTE2: When evolve is called from an island, the process is forked and transferred to another python or ipython
instance. As a consequence, when writing your _obj_fun_impl you cannot use stuff like matplotlib to 
make interactive plots and alike. If you need, during development, to have this kind of support,
use the algorithm evolve method, for example

.. code-block:: python

   from PyGMO import *
   from  my_module import my_problem

   prob = my_problem(dim=10)
   algo = algorithm.bee_colony(gen=100)
   isl = island(algo,prob,20)
   pop = island.population
   pop = algo.evolve(pop)
   print pop.champion.f

NOTE3: If performance is your goal, you should implement your problem in C++, and then expose it into python.