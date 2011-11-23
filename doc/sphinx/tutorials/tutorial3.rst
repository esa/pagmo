.. _tutorial3:

================================================================
Tutorial 3: Coding your own algorithm in Python
================================================================

In this tutorial we will learn how to code a simple optimization algorithm.
We will write the algorithm so that it manage multi-objective, mixed_int, constrained optimization
as this will allow us to explain all the basic PyGMO functioning. Clearly our algorithm will not
be very good ... a random search always useful to benchmark against :)

In a nutshell ... we will write a class deriving from :class:`PyGMO.algorithm.base`
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

	#This is the 'juice' of the algorithm, the method where the actual optimzation is coded. 
	def evolve(self,pop):
		#If the population is empty (i.e. no individuals) nothing happens
		if len(pop) == 0:
			return pop
			
		#Here we rename some variables, in particular the problem
		prob = pop.problem
		#Its dimensions (total and continuous)
		dim, cont_dim = prob.dimension, prob.dimension - prob.i_dimension
		#And the lower/upper bounds for the chromosome
		lb, ub = prob.lb, prob.ub
		import random

		#The algorithm now starts manipulating the population
		for _ in range(self.__iter):
			#We create a random vector within the bounds ... first the continuous part
			tmp_cont = [random.uniform(lb[i],ub[i]) for i in range(cont_dim)]
			#then the integer part
			tmp_int = [float(random.randint(lb[i],ub[i])) for i in range(cont_dim,dim)]
			#and we assemble them into one decision vector
			tmp_x = tmp_cont + tmp_int
			#we compute the objective function of our mutated vector
			tmp_f = prob.objfun(tmp_x)
			#and the value of the constraints
			tmp_c = prob.compute_constraints(tmp_x)
			#we extract the current worst population individual
			worst_idx = pop.get_worst_idx()
			worst = pop[worst_idx]
			#and subsitute it with the mutated if this is actually better
			if prob.compare_fc(tmp_f,tmp_c,worst.cur_f,worst.cur_c):
				pop.set_x(worst_idx,tmp_x)
		#at the end of it all we return the 'evolved' population
		return pop

	def get_name(self):
		return "Monte Carlo (Python)"

	def human_readable_extra(self):
		return "n_iter=" + str(self.__n_iter)

The above code contains a lot of interesting points worth to be discussed. So, we start

* In PyGMO the decision vector (chromosome) is represented as an n-tuple. Its dimension and structure depends
  on the problem. Its dimension will be problem.dimension, the first prob.dimension - prob.i_dimension components will
  be continuous, the remaining problem.i_dimension will instead be integers.
* The method prob.compute_constraints is virtual, its default implementation returns an empty tuple. This allow
  our algorithm to 'work' on constrained as well unconstrained problems.
* The method prob.compare_fc is virtual. Its default implementation counts the number of satisfied constraints, then
  the number of dominated individuals. This allow our algorithm to work for single, as well as for multi-objective problems
