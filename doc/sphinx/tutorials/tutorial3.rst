================================================================
Tutorial 3: Coding your own algorithm in Python
================================================================

In this tutorial we will learn how to code a simple optimization algorithm.
We will write the algorithm so that it manage multi-objective, mixed_int, constrained optimization
as this will allow us to explain all the basic PyGMO functioning. Clearly our algorithm will not
be very good ... a random search always useful to benchmark against :)

In a nutshell ... we will write a class deriving from algorithm.base
and reimplement some of its 'virtual' methods, the main one being evolve!!!

.. code-block:: python

   class my_algorithm(base):
	"""
	Monte-Carlo (random sampling) algorithm implemented purely in Python.
	"""

	def __init__(self,iter = 10):
		"""
		Constructs a Monte-Carlo (random sampling) algorithm
		
		USAGE: algorithm.my_algorithm(iter = 10)
		
		NOTE: At the end of each iteration, the randomly generated 
			point substitutes the worst individual in the population if better
		
		* iter: number of random samples
		"""
		#We start calling the base constructor
		super(my_algorithm,self).__init__()
		#We then define the algorithm 'private' data members
		self.__iter = iter

	def evolve(self,pop):
		if len(pop) == 0:
			return pop
		prob = pop.problem
		dim, cont_dim = prob.dimension, prob.dimension - prob.i_dimension
		lb, ub = prob.lb, prob.ub
		import random
		for _ in range(self.__n_iter):
			tmp_cont = [random.uniform(lb[i],ub[i]) for i in range(cont_dim)]
			tmp_int = [float(random.randint(lb[i],ub[i])) for i in range(cont_dim,dim)]
			tmp_x = tmp_cont + tmp_int
			tmp_f = prob.objfun(tmp_x)
			tmp_c = prob.compute_constraints(tmp_x)
			worst_idx = pop.get_worst_idx()
			worst = pop[worst_idx]
			if prob.compare_fc(tmp_f,tmp_c,worst.cur_f,worst.cur_c):
				pop.set_x(worst_idx,tmp_x)
		return pop

	def get_name(self):
		return "Monte Carlo (Python)"

	def human_readable_extra(self):
		return "n_iter=" + str(self.__n_iter)