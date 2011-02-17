# -*- coding: utf-8 -*-
from _algorithm import *

# Raw C++ base class.
_base = _algorithm._base

class base(_algorithm._base):
	def __init__(self):
		_algorithm._base.__init__(self)
	def get_name(self):
		return str(type(self))
	def __get_deepcopy__(self):
		from copy import deepcopy
		return deepcopy(self)

class py_test(base):
	"""
	Simple Monte Carlo algorithm.
	"""
	def __init__(self,n_iter = 10):
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

# Helper class to ease the inclusion of scipy.optimize solvers.
class _scipy_base(base):
	def __init__(self,solver_name,constrained,maxiter,tol):
		base.__init__(self)
		try:
			exec('from scipy.optimize import %s as solver' % solver_name)
			from numpy import concatenate, array
		except ImportError:
			raise ImportError('The necessary SciPy and/or NumPy classes/functions could not be imported')
		self.solver = solver
		self.constrained = constrained
		self.maxiter = maxiter
		self.tol = tol
	# Check if problem is compatible with the algorithm.
	def _problem_checks(self,prob):
		if prob.f_dimension > 1:
			raise ValueError("this algorithm does not support multi-objective optimisation")
		if prob.dimension == prob.i_dimension:
			raise ValueError("the provided problem has no continuous part")
		if not self.constrained and prob.c_dimension:
			raise ValueError("this algorithm does not support constrained optimisation")
	# Check that the algorithm did not go out of bounds, and, in such a case, correct the chromosome.
	def _check_new_chromosome(self,new_chromosome,prob):
		for i in range(0,len(new_chromosome)):
			if new_chromosome[i] < prob.lb[i]:
				new_chromosome[i] = prob.lb[i]
			if new_chromosome[i] > prob.ub[i]:
				new_chromosome[i] = prob.ub[i]
		return new_chromosome
	def _starting_params(self,pop):
		from numpy import array
		# Number of equality constraints.
		n_ec = pop.problem.c_dimension - pop.problem.ic_dimension
		# Extract the continuous part of the first individual's current chromosome.
		x0 = array(pop[pop.get_best_idx()].cur_x[0:pop.problem.dimension - pop.problem.i_dimension],dtype=float)
		# Combinatorial part of the chromosome (which will not be optimised).
		x0_comb = array(pop[pop.get_best_idx()].cur_x[pop.problem.dimension - pop.problem.i_dimension:],dtype=float)
		return n_ec,x0,x0_comb
	def _human_readable_extra(self):
		return "max_iter = " + str(self.maxiter) + ", tol = " + str(self.tol) + "\n"

class scipy_fmin(_scipy_base):
	"""
	Wrapper around SciPy's fmin optimiser.
	"""
	def __init__(self,maxiter = 100,tol = 1E-5,verbose = False):
		_scipy_base.__init__(self,'fmin',False,maxiter,tol)
		self.verbose = verbose
	def evolve(self,pop):
		from numpy import concatenate
		prob = pop.problem
		self._problem_checks(prob)
		if len(pop) == 0:
			return pop
		_, x0, x0_comb = self._starting_params(pop)
		retval = self.solver(lambda x: prob.objfun(concatenate((x,x0_comb)))[0],x0,disp = int(self.verbose),ftol = self.tol,maxiter = self.maxiter)
		new_chromosome = list(retval) + list(x0_comb)
		pop.set_x(pop.get_best_idx(),self._check_new_chromosome(new_chromosome,prob))
		return pop

class scipy_l_bfgs_b(_scipy_base):
	"""
	Wrapper around SciPy's l_bfgs_b optimiser.
	"""
	def __init__(self,maxiter = 100,tol = 1E-5,verbose = False):
		_scipy_base.__init__(self,'fmin_l_bfgs_b',False,maxiter,tol)
		self.verbose = verbose
	def evolve(self,pop):
		from numpy import concatenate, array
		prob = pop.problem
		self._problem_checks(prob)
		if len(pop) == 0:
			return pop
		_, x0, x0_comb = self._starting_params(pop)
		# Extract the a list of tuples representing the bounds.
		prob_bounds = [(prob.lb[i],prob.ub[i]) for i in range(0,prob.dimension - prob.i_dimension)]
		if self.verbose:
			iprn = 0
		else:
			iprn = -1
		retval = self.solver(lambda x: array(prob.objfun(concatenate((x,x0_comb))),dtype=float),x0,bounds = prob_bounds,approx_grad = True, iprint = iprn,pgtol = self.tol, maxfun = self.maxiter)
		new_chromosome = list(retval[0]) + list(x0_comb)
		pop.set_x(pop.get_best_idx(),self._check_new_chromosome(new_chromosome,prob))
		return pop

class scipy_slsqp(_scipy_base):
	"""
	Wrapper around SciPy's slsqp optimiser.
	"""
	def __init__(self,maxiter = 100,tol = 1E-6,verbose = False):
		_scipy_base.__init__(self,'fmin_slsqp',True,maxiter,tol)
		self.verbose = verbose
	def get_name(self):
		return 'fmin_slsqp'
	def evolve(self,pop):
		from numpy import concatenate, array
		prob = pop.problem
		self._problem_checks(prob)
		# If population is empty, just return input population.
		if len(pop) == 0:
			return pop
		# Get starting params.
		n_ec, x0, x0_comb = self._starting_params(pop)
		# Extract the a list of tuples representing the bounds.
		prob_bounds = [(prob.lb[i],prob.ub[i]) for i in range(0,prob.dimension - prob.i_dimension)]
		if self.verbose:
			iprn = 1
		else:
			iprn = 0
		# Run the optimisation.
		retval = self.solver(lambda x: prob.objfun(concatenate((x, x0_comb)))[0],x0,f_eqcons = lambda x: array(prob.compute_constraints(concatenate((x, x0_comb)))[0:n_ec],dtype=float),
			f_ieqcons = lambda x: array(prob.compute_constraints(concatenate((x, x0_comb)))[n_ec:],dtype=float) * -1,bounds = prob_bounds,iprint = iprn,
			iter = self.maxiter,acc = self.tol)
		# Set the individual's chromosome in the population and return. Conserve the integer part from the
		# original individual.
		new_chromosome = list(retval) + list(x0_comb)
		pop.set_x(pop.get_best_idx(),self._check_new_chromosome(new_chromosome,prob))
		return pop

class scipy_tnc(_scipy_base):
	"""
	Wrapper around SciPy's tnc optimiser.
	"""
	def __init__(self,maxiter = 100,tol = 1E-6,verbose = False):
		_scipy_base.__init__(self,'fmin_tnc',False,maxiter,tol)
		self.verbose = verbose
	def evolve(self,pop):
		from numpy import concatenate, array
		prob = pop.problem
		self._problem_checks(prob)
		if len(pop) == 0:
			return pop
		_, x0, x0_comb = self._starting_params(pop)
		# Extract the a list of tuples representing the bounds.
		prob_bounds = [(prob.lb[i],prob.ub[i]) for i in range(0,prob.dimension - prob.i_dimension)]
		if self.verbose:
			msg = 15
		else:
			msg = 0
		retval = self.solver(lambda x: array(prob.objfun(concatenate((x,x0_comb))),dtype=float),x0,bounds = prob_bounds,approx_grad = True, messages = msg,
			pgtol = self.tol, maxfun = self.maxiter)
		new_chromosome = list(retval[0]) + list(x0_comb)
		pop.set_x(pop.get_best_idx(),self._check_new_chromosome(new_chromosome,prob))
		return pop

class scipy_cobyla(_scipy_base):
	"""
	Wrapper around SciPy's cobyla optimiser.
	"""
	def __init__(self,maxiter = 100,tol = 1E-5,verbose = False):
		_scipy_base.__init__(self,'fmin_cobyla',True,maxiter,tol)
		self.verbose = verbose
	def evolve(self,pop):
		from numpy import concatenate, array
		prob = pop.problem
		self._problem_checks(prob)
		if len(pop) == 0:
			return pop
		nec, x0, x0_comb = self._starting_params(pop)
		# We need to build a vector of functions for the constraints. COBYLA uses inequality constraints with >= 0.
		ub = prob.ub
		lb = prob.lb
		f_cons = []
		for i in range(0,nec):
			# Each eq. constraint is converted into two ineq. constraints with different signs,
			f_cons.append(lambda x, i = i: prob.compute_constraints(concatenate((x, x0_comb)))[i])
			f_cons.append(lambda x, i = i: -prob.compute_constraints(concatenate((x, x0_comb)))[i])
		for i in range(nec,prob.c_dimension):
			# Ineq. constraints.
			f_cons.append(lambda x, i = i: -prob.compute_constraints(concatenate((x, x0_comb)))[i])
		for i in range(0,prob.dimension - prob.i_dimension):
			# Box bounds implemented as inequality constraints.
			f_cons.append(lambda x, i = i: x[i] - lb[i])
			f_cons.append(lambda x, i = i: ub[i] - x[i])
		if self.verbose:
			iprn = 1
		else:
			iprn = 0
		retval = self.solver(lambda x: array(prob.objfun(concatenate((x,x0_comb))),dtype=float),x0,cons = f_cons,iprint = iprn, maxfun = self.maxiter, rhoend = self.tol)
		new_chromosome = list(retval) + list(x0_comb)
		pop.set_x(pop.get_best_idx(),self._check_new_chromosome(new_chromosome,prob))
		return pop

class scipy_anneal(_scipy_base):
	"""
	Wrapper around SciPy's anneal optimiser.
	"""
	def __init__(self,maxiter = 100,tol = 1E-5,verbose = False):
		_scipy_base.__init__(self,'anneal',False,maxiter,tol)
		self.verbose = verbose
	def evolve(self,pop):
		from numpy import concatenate, array
		prob = pop.problem
		self._problem_checks(prob)
		# If population is empty, just return input population.
		if len(pop) == 0:
			return pop
		# Get starting params.
		n_ec, x0, x0_comb = self._starting_params(pop)
		# Run the optimisation.
		retval = self.solver(lambda x: prob.objfun(concatenate((x, x0_comb)))[0],x0,lower = array(prob.lb,dtype=float),upper = array(prob.ub,dtype=float),
			full_output = int(self.verbose), maxiter = self.maxiter, feps = self.tol)
		# Set the individual's chromosome in the population and return. Conserve the integer part from the
		# original individual.
		new_chromosome = list(retval[0]) + list(x0_comb)
		pop.set_x(0,self._check_new_chromosome(new_chromosome,prob))
		return pop

def _get_algorithm_list():
	from PyGMO import algorithm
	# Try importing SciPy and NumPy.
	try:
		import scipy, numpy
		algorithm_list = [algorithm.__dict__[n] for n in filter(lambda n: not n.startswith('_') and not n == 'base',dir(algorithm))]
	except ImportError as e:
		algorithm_list = [algorithm.__dict__[n] for n in filter(lambda n: not n.startswith('_') and not n == 'base' and not n.startswith('scipy'),dir(algorithm))]
	return algorithm_list
