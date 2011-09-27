from base import base

class py_example(base):
	"""
	Monte-Carlo Algorithm (Python)
	"""
	def __init__(self,n_iter = 10):
		"""
		Constructs a Monte Carlo Algorithm (used as an example of algorithm building directly in Python)
		
		USAGE: algorithm.py_example(iter = 10000)
		
		NOTE: At the end of each iteration, the randomly generated 
			point substitutes the worst in the population if better
		
		* iter: number of Monte Carlo runs
		"""
		base.__init__(self)
		self.__n_iter = n_iter
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