# -*- coding: utf-8 -*-
from _algorithm import *
from base import base
from example import py_example
from cross_entropy import py_cross_entropy

_base = _algorithm._base

# Renaming and placing the enums
_algorithm.sga.crossover = _algorithm._crossover_type
_algorithm.sga.selection = _algorithm._selection_type
_algorithm.sga.mutation = _algorithm._mutation_type

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
		

def _de_ctor(self,**kwargs):
	"""
	Constructs a Differential Evolution algorithm:
	
	USAGE: algorithm.de(gen=1,f=0.8,cr=0.9,variant=2)
	
	* gen: number of generations
	* f: weighting factor in [0,1]
	* cr: crossover in [0,1]
	* variant: algoritmic variant to use (one of [1 .. 10])
		1. DE/best/1/exp
		2. DE/rand/1/exp
		3. DE/rand-to-best/1/exp
		4. DE/best/2/exp
		5. DE/rand/2/exp
		6. DE/best/1/bin
		7. DE/rand/1/bin
		8. DE/rand-to-best/1/bin
		9. DE/best/2/bin
		10. DE/rand/2/bin
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(kwargs.pop('gen', 1))
	arg_list.append(kwargs.pop('f', 0.8))
	arg_list.append(kwargs.pop('cr', 0.9))
	arg_list.append(kwargs.pop('variant', 2))
	self._orig_init(*arg_list)
de._orig_init = de.__init__
de.__init__ = _de_ctor

def _pso_ctor(self,**kwargs):
	"""
	Constructs a Particle Swarm Optimization (steady-state). The position update is applied
	immediately after the velocity update
	
	USAGE: algorithm.pso(gen=1, omega = 0.7298, eta1 = 2.05, eta2 = 2.05, vcoeff = 0.5, variant = 5, neighb_type = 2, neighb_param = 4])

	* gen: number of generations
	* omega: constriction factor (or particle inertia weight) in [0,1]
	* eta1: Cognitive component in [0,4]
	* eta2: Social component in [0,4]
	* vcoeff: Maximum velocity coefficient (w.r.t. the box-bounds width) in [0,1]
	* variant: algoritmic variant to use (one of  [1 .. 6])
		1. PSO canonical (with inertia weight)
		2. PSO canonical (with inertia weight
			and equal random weights of social and cognitive components)
		3. PSO variant (with inertia weight
			same random number for all components.)
		4. PSO variant (with inertia weight
			same random number for all components
			and equal weights of social and cognitive components)
		5. PSO canonical (with constriction factor)
		6. Fully Informed Particle Swarm (FIPS)
	* neighb_type: defines the particle neighbourhood (used for the social component)
		1. gbest neighbourhood topology (fully connected)
		2. lbest neighbourhood topology (ring)
		3. Von-Neumann neighbourhood topology (square lattice)
		4. Randomly-varying neighbourhood topology (not yet implemented)
	* neighb_param: in the lbest topology defines how many links 
		are there in the ring (half to the right and half to the left)
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(kwargs.pop('gen', 1))
	arg_list.append(kwargs.pop('omega', 0.7298))
	arg_list.append(kwargs.pop('eta1', 2.05))
	arg_list.append(kwargs.pop('eta2', 2.05))
	arg_list.append(kwargs.pop('vcoeff', 0.5))
	arg_list.append(kwargs.pop('variant', 5))
	arg_list.append(kwargs.pop('neighb_type', 2))
	arg_list.append(kwargs.pop('neighb_param', 4))	
	self._orig_init(*arg_list)
pso._orig_init = pso.__init__
pso.__init__ = _pso_ctor

def _pso_gen_ctor(self,**kwargs):
	"""
	Constructs a Particle Swarm Optimization (generational). The position update is applied
	immediately after the velocity update
	
	USAGE: algorithm.pso_gen(gen=1, omega = 0.7298, eta1 = 2.05, eta2 = 2.05, vcoeff = 0.5, variant = 5, neighb_type = 2, neighb_param = 4])

	* gen: number of generations
	* omega: constriction factor (or particle inertia weight) in [0,1]
	* eta1: Cognitive component in [0,4]
	* eta2: Social component in [0,4]
	* vcoeff: Maximum velocity coefficient (w.r.t. the box-bounds width) in [0,1]
	* variant: algoritmic variant to use (one of  [1 .. 6])
		1. PSO canonical (with inertia weight)
		2. PSO canonical (with inertia weight
			and equal random weights of social and cognitive components)
		3. PSO variant (with inertia weight
			same random number for all components.)
		4. PSO variant (with inertia weight
			same random number for all components
			and equal weights of social and cognitive components)
		5. PSO canonical (with constriction factor)
		6. Fully Informed Particle Swarm (FIPS)
	* neighb_type: defines the particle neighbourhood (used for the social component)
		1. gbest neighbourhood topology (fully connected)
		2. lbest neighbourhood topology (ring)
		3. Von-Neumann neighbourhood topology (square lattice)
		4. Randomly-varying neighbourhood topology (not yet implemented)
	* neighb_param: in the lbest topology defines how many links 
		are there in the ring (half to the right and half to the left)
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(kwargs.pop('gen', 1))
	arg_list.append(kwargs.pop('omega', 0.7298))
	arg_list.append(kwargs.pop('eta1', 2.05))
	arg_list.append(kwargs.pop('eta2', 2.05))
	arg_list.append(kwargs.pop('vcoeff', 0.5))
	arg_list.append(kwargs.pop('variant', 5))
	arg_list.append(kwargs.pop('neighb_type', 2))
	arg_list.append(kwargs.pop('neighb_param', 4))	
	self._orig_init(*arg_list)
pso_gen._orig_init = pso_gen.__init__
pso_gen.__init__ = _pso_gen_ctor

def _sga_ctor(self,**kwargs):
	"""
	Constructs a Simple Genetic Algorithm (generational)
	
	USAGE: algorithm.sga(gen=1, cr=.95, m=.02, elitism=1, mutation=GAUSSIAN, width = 0.1, selection=ROULETTE, crossover=EXPONENTIAL)
  
	* gen: number of generations
	* cr: crossover factor in [0,1]
	* m: mutation probability (for each component) [0,1]
	* elitism: number of generation after which the best is reinserted
	* mutation: mutation type (one of [RANDOM, GAUSSIAN])
	* width: the mutation width (in case of a GAUSSIAN bell
		this is the std normalized with the width)
	* selection: selection startegy (one of [ROULETTE, BEST20])
	* crossover: crossover strategy (one of [BINOMIAL, EXPONENTIAL])
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(kwargs.pop('gen', 1))
	arg_list.append(kwargs.pop('cr', 0.95))
	arg_list.append(kwargs.pop('m', 0.02))
	arg_list.append(kwargs.pop('elitism', 1))
	arg_list.append(kwargs.pop('mutation', sga.mutation.GAUSSIAN))
	arg_list.append(kwargs.pop('width', 0.1))
	arg_list.append(kwargs.pop('selection', sga.selection.ROULETTE))
	arg_list.append(kwargs.pop('crossover', sga.crossover.EXPONENTIAL))	
	self._orig_init(*arg_list)
sga._orig_init = sga.__init__
sga.__init__ = _sga_ctor

def _sa_corana_ctor(self,**kwargs):
	"""
	Constructs Corana's Simulated Annealing
	
	USAGE: algorithm.sa_corana(iter = 10000, Ts = 10, Tf = .1, steps = 1, bin_size = 20, range = 1)
	
	NOTE: as this version of simulated annealing loops through the chromosome, the iter number needs to be selected
	large enough to allow the temperature schedule to actuallt make sense. For example if your problem has D dimensions
	then in order to have at least N temperature adjustments (from Ts to Tf) one should select iter = D * N * steps * bin_size.

	* iter: number of total iterations
	* Ts: starting temperature
	* Tf: final temperature ( > Ts)
	* steps: number of steps adjustments
	* bin_size: size of the bin used to evaluate the step adjustment
	* range: initial size of the neighbourhood (in [0,1])
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(kwargs.pop('iter', 10000))
	arg_list.append(kwargs.pop('Ts', 10))
	arg_list.append(kwargs.pop('Tf', 0.1))
	arg_list.append(kwargs.pop('steps', 1))
	arg_list.append(kwargs.pop('bin_size', 20))
	arg_list.append(kwargs.pop('range', 1))	
	self._orig_init(*arg_list)
sa_corana._orig_init = sa_corana.__init__
sa_corana.__init__ = _sa_corana_ctor

def _bee_colony_ctor(self,**kwargs):
	"""
	Constructs an Artificial Bee Colony Algorithm
	
	USAGE: algorithm.bee_colony(gen = 1, limit = 20)
	
	* gen: number of 'generations' (each generation 2*NP function evaluations
		are made where NP is the population size)
	* limit: number of tries after which a source of food is dropped if not improved
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(kwargs.pop('gen', 1))
	arg_list.append(kwargs.pop('limit', 20))
	self._orig_init(*arg_list)
bee_colony._orig_init = bee_colony.__init__
bee_colony.__init__ = _bee_colony_ctor

def _firefly_ctor(self,**kwargs):
	"""
	Constructs a Firefly Algorithm
	
	USAGE: algorithm.firefly(gen = 1, alpha = 0.01, beta = 1.0, gamma = 0.8)
	
	* gen: number of 'generations' 
	* alpha: width of the random vector (in [0,1])
	* beta: maximum attractiveness (in [0,1])
	* gamma: absorption coefficient (in [0,1])
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(kwargs.pop('gen', 1))
	arg_list.append(kwargs.pop('alpha', 20))
	arg_list.append(kwargs.pop('beta', 20))
	arg_list.append(kwargs.pop('gamma', 20))
	self._orig_init(*arg_list)
firefly._orig_init = firefly.__init__
firefly.__init__ = _firefly_ctor

def _ms_ctor(self,**kwargs):
	"""
	Constructs a Multistart Algorithm
	
	USAGE: algorithm.ms(algorithm = algorithm.de(), iter = 1)
	
	NOTE: starting from pop1, at each iteration a random pop2 is evolved
	with the selected algorithm and its final best replaces the worst of pop1

	* algorithm: PyGMO algorithm to be multistarted
	* iter: number of multistarts

	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(kwargs.pop('algorithm', de()))
	arg_list.append(kwargs.pop('iter', 1))
	self._orig_init(*arg_list)
ms._orig_init = ms.__init__
ms.__init__ = _ms_ctor

def _mbh_ctor(self,**kwargs):
	"""
	Constructs a Monotonic Basin Hopping Algorithm (generalized to accept any algorithm)
	
	USAGE: algorithm.mbh(algorithm = algorithm.cs(), stop = 5, perturb = 5e-2);
	
	NOTE: Starting from pop, algorithm is applied to the perturbed pop returning pop2. If pop2 is better than
	pop then pop=pop2 and a counter is reset to zero. If pop2 is not better the counter is incremented. If
	the counter is larger than stop, optimization is terminated

	* algorithm: 'local' optimiser
	* stop: number of no improvements before halting the optimization
	* perturb: non-dimentional perturbation width (can be a list, in which case
		it has to have the same dimension of the problem mbh will be applied to)
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(kwargs.pop('algorithm', de()))
	arg_list.append(kwargs.pop('stop', 1))
	arg_list.append(kwargs.pop('perturb', 5e-2))
	self._orig_init(*arg_list)
mbh._orig_init = mbh.__init__
mbh.__init__ = _mbh_ctor

def _cs_ctor(self,**kwargs):
	"""
	Constructs a Compass Search Algorithm 
	
	USAGE: algorithm.cs(max_eval = 1, stop_range = 0.01, start_range = 0.1, reduction_coeff = 0.5);
	

	* max_iter: maximum number of function evaluations
	* stop_range: when the range is reduced to a value smaller than stop_range cs stops
	* start_range: starting range (non-dimensional wrt ub-lb)
	* reduction_coeff: the range is multiplied by reduction_coeff whenever no improvment is made
	                   across one chromosome
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(kwargs.pop('max_eval', 10000))
	arg_list.append(kwargs.pop('stop_range', 0.01))
	arg_list.append(kwargs.pop('start_range', 0.1))
	arg_list.append(kwargs.pop('reduction_coeff', 0.5))
	self._orig_init(*arg_list)
cs._orig_init = cs.__init__
cs.__init__ = _cs_ctor

def _ihs_ctor(self,**kwargs):
	"""
	Constructs an Improved Harmony Search Algorithm
	
	USAGE: algorithm.ihs(gen = 1, hmcr = 0.85, par_min = 0.35, par_max = 0.99, bw_min = 1E-5, bw_max = 1);
	
	* gen: number of generations (improvisations)
	* hmcr: rate of choosing from memory (in ]0,1[)
	* par_min: minimum pitch adjustment rate (in ]0,1[)
	* par_max: maximum pitch adjustment rate (in ]0,1[, > par_min) 
	* bw_min: minimum distance bandwidth 
	* bw_max: maximum distance bandwidth (> bw_min)
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(kwargs.pop('gen', 1))
	arg_list.append(kwargs.pop('hmcr', 0.85))
	arg_list.append(kwargs.pop('par_min', 0.35))
	arg_list.append(kwargs.pop('par_max', 0.99))
	arg_list.append(kwargs.pop('bw_min', 1e-5))
	arg_list.append(kwargs.pop('bw_max', 1))
	self._orig_init(*arg_list)
ihs._orig_init = ihs.__init__
ihs.__init__ = _ihs_ctor

def _monte_carlo_ctor(self,**kwargs):
	"""
	Constructs a Monte Carlo Algorithm
	
	USAGE: algorithm.monte_carlo(iter = 10000)
	
	NOTE: At the end of each iteration, the randomly generated 
		point substitutes the worst in the population if better
	
	* iter: number of Monte Carlo runs
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(kwargs.pop('iter', 10000))
	self._orig_init(*arg_list)
monte_carlo._orig_init = monte_carlo.__init__
monte_carlo.__init__ = _monte_carlo_ctor

def _get_algorithm_list():
	from PyGMO import algorithm
	# Try importing SciPy and NumPy.
	try:
		import scipy, numpy
		algorithm_list = [algorithm.__dict__[n] for n in filter(lambda n: not n.startswith('_') and not n == 'base',dir(algorithm))]
	except ImportError as e:
		algorithm_list = [algorithm.__dict__[n] for n in filter(lambda n: not n.startswith('_') and not n == 'base' and not n.startswith('scipy'),dir(algorithm))]
	return algorithm_list
