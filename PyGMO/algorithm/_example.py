from _base import base

class py_example(base):
	"""
	Monte-Carlo (random sampling) algorithm implemented purely in Python.
	"""

	def __init__(self,iter = 10):
		"""
		Constructs a Monte-Carlo (random sampling) algorithm
		
		USAGE: algorithm.py_example(iter = 10)
		
		NOTE: At the end of each iteration, the randomly generated 
			point substitutes the worst individual in the population if better
		
		* iter: number of random samples
		"""
		#We start calling the base constructor
		super(py_example,self).__init__()
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
		return "iter=" + str(self.__iter)

