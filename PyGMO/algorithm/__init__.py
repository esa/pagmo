# -*- coding: iso-8859-1 -*-
from _algorithm import *

# Helper class to ease the inclusion of scipy.optimize solvers.
class _scipy_optimize_algorithm(object):
	def __init__(self,solver_name,constrained = False):
		try:
			exec('from scipy.optimize import %s as solver' % solver_name)
			from numpy import concatenate, array
		except ImportError:
			raise ImportError('The necessary SciPy and/or NumPy classes/functions could not be imported')
		self.solver = solver
		self.constrained = constrained
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
		x0 = array(pop[0].cur_x[0:pop.problem.dimension - pop.problem.i_dimension],dtype=float)
		# Combinatorial part of the chromosome (which will not be optimised).
		x0_comb = array(pop[0].cur_x[pop.problem.dimension - pop.problem.i_dimension:],dtype=float)
		return n_ec,x0,x0_comb

class scipy_fmin(base,_scipy_optimize_algorithm):
	"""
	Wrapper around SciPy's fmin optimiser.
	"""
	def __init__(self,verbose = False):
		base.__init__(self)
		_scipy_optimize_algorithm.__init__(self,'fmin',constrained = False)
		self.verbose = verbose
	def __copy__(self):
		return scipy_fmin(verbose = self.verbose)
	def evolve(self,pop):
		from numpy import concatenate
		prob = pop.problem
		self._problem_checks(prob)
		if len(pop) == 0:
			return pop
		_, x0, x0_comb = self._starting_params(pop)
		retval = self.solver(lambda x: prob.objfun(concatenate((x,x0_comb)))[0],x0,disp = int(self.verbose))
		new_chromosome = list(retval) + list(x0_comb)
		pop.set_x(0,self._check_new_chromosome(new_chromosome,prob))
		return pop

class scipy_l_bfgs_b(base,_scipy_optimize_algorithm):
	"""
	Wrapper around SciPy's l_bfgs_b optimiser.
	"""
	def __init__(self,verbose = False):
		base.__init__(self)
		_scipy_optimize_algorithm.__init__(self,'fmin_l_bfgs_b',constrained = False)
		self.verbose = verbose
	def __copy__(self):
		return scipy_l_bfgs_b(self.verbose)
	def evolve(self,pop):
		from numpy import concatenate, array
		prob = pop.problem
		self._problem_checks(prob)
		if len(pop) == 0:
			return pop
		_, x0, x0_comb = self._starting_params(pop)
		# Extract the a list of tuples representing the bounds.
		prob_bounds = [(prob.lb[i],prob.ub[i]) for i in range(0,prob.dimension)]
		if self.verbose:
			iprn = 0
		else:
			iprn = -1
		retval = self.solver(lambda x: array(prob.objfun(concatenate((x,x0_comb))),dtype=float),x0,bounds = prob_bounds,approx_grad = True, iprint = iprn)
		new_chromosome = list(retval[0]) + list(x0_comb)
		pop.set_x(0,self._check_new_chromosome(new_chromosome,prob))
		return pop

class scipy_slsqp(base,_scipy_optimize_algorithm):
	"""
	Wrapper around SciPy's slsqp optimiser.
	"""
	def __init__(self,verbose = False):
		base.__init__(self)
		_scipy_optimize_algorithm.__init__(self,'fmin_slsqp',constrained = True)
		self.verbose = verbose
	def __copy__(self):
		return scipy_slsqp(self.verbose)
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
		prob_bounds = [(prob.lb[i],prob.ub[i]) for i in range(0,prob.dimension)]
		if self.verbose:
			iprn = 1
		else:
			iprn = 0
		# Run the optimisation.
		retval = self.solver(lambda x: prob.objfun(concatenate((x, x0_comb)))[0],x0,f_eqcons = lambda x: array(prob.compute_constraints(concatenate((x, x0_comb)))[0:n_ec],dtype=float),
			f_ieqcons = lambda x: array(prob.compute_constraints(concatenate((x, x0_comb)))[n_ec:],dtype=float) * -1,bounds = prob_bounds,iprint = iprn)
		# Set the individual's chromosome in the population and return. Conserve the integer part from the
		# original individual.
		new_chromosome = list(retval) + list(x0_comb)
		pop.set_x(0,self._check_new_chromosome(new_chromosome,prob))
		return pop
