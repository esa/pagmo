.. _adding_a_new_algorithm:

================================================================
Adding a new algorithm
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
			#which we push back in the population
			pop.push_back(tmp_x)
			#to then remove the worst individual
			pop.erase(pop.get_worst_idx())
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
* When a chromosome is pushed back in a population, the domination count and the domination list (data members)
  are automatically updated
* The get_worst_idx method sort the population with respect to the domination count, then domination
  list size (in inverse order)
