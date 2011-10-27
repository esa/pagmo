# -*- coding: utf-8 -*-
from _algorithm import *
from _algorithm import _base
from _base import base
from _example import py_example
from _cross_entropy import py_cross_entropy
from _scipy_algos import *

_base = _algorithm._base

# Renaming and placing the enums
_algorithm.sga.crossover = _algorithm._crossover_type
_algorithm.sga.selection = _algorithm._selection_type
_algorithm.sga.mutation = _algorithm._mutation_type

#Creating the list of algorithms
def _get_algorithm_list():
	import _algorithm as algorithm
	# Try importing SciPy and NumPy.
	try:
		import scipy, numpy
		algorithm_list = [algorithm.__dict__[n] for n in filter(lambda n: not n.startswith('_') and not n == 'base',dir(algorithm))]
		from PyGMO import __version__
		__version__ = __version__ + "GTOP " + "GSL " + "SciPy "
	except ImportError as e:
		algorithm_list = [algorithm.__dict__[n] for n in filter(lambda n: not n.startswith('_') and not n == 'base' and not n.startswith('scipy'),dir(algorithm))]
	return algorithm_list

# Redefining the constructors of all algorithms to obtain good documentation and to allow kwargs
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
		4. Randomly-varying neighbourhood topology
	* neighb_param: if the lbest topology is selected, it represents each particle's indegree
		(also outdegree) in the swarm topology. Particles have neighbours up
		to a radius of k = neighb_param / 2 in the ring. If the Randomly-varying neighbourhood topology
		is selected, neighb_param represents each particle's maximum outdegree in the swarm topology.
		The minimum outdegree is 1 (the particle always connects back to itself).
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
		4. Randomly-varying neighbourhood topology
	* neighb_param: if the lbest topology is selected, it represents each particle's indegree
		(also outdegree) in the swarm topology. Particles have neighbours up
		to a radius of k = neighb_param / 2 in the ring. If the Randomly-varying neighbourhood topology
		is selected, neighb_param represents each particle's maximum outdegree in the swarm topology.
		The minimum outdegree is 1 (the particle always connects back to itself).
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

#def _firefly_ctor(self,**kwargs):
#	"""
#	Constructs a Firefly Algorithm
#	
#	USAGE: algorithm.firefly(gen = 1, alpha = 0.01, beta = 1.0, gamma = 0.8)
#	
#	* gen: number of 'generations' 
#	* alpha: width of the random vector (in [0,1])
#	* beta: maximum attractiveness (in [0,1])
#	* gamma: absorption coefficient (in [0,1])
#	"""
#	# We set the defaults or the kwargs
#	arg_list=[]
#	arg_list.append(kwargs.pop('gen', 1))
#	arg_list.append(kwargs.pop('alpha', 20))
#	arg_list.append(kwargs.pop('beta', 20))
#	arg_list.append(kwargs.pop('gamma', 20))
#	self._orig_init(*arg_list)
#firefly._orig_init = firefly.__init__
#firefly.__init__ = _firefly_ctor

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
	

	* max_eval: maximum number of function evaluations
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

#NLOPT algorithms (only if PyGMO has been compiled woth nlopt option activated)
if "nlopt" in str(_get_algorithm_list()):
	def _nlopt_bobyqa_ctor(self,**kwargs):
		"""
		Constructs a BOBYQA algorithm (Bound Optimization BY Quadratic Approximation) (NLOPT)
	
		USAGE: algorithm.nlopt_bobyqa(max_iter = 100, ftol = 1e-6, xtol = 1e-6);
	
		* max_iter: stop-criteria (number of iterations)
		* ftol: stop-criteria (absolute on the obj-fun)
		* xtol: stop-criteria (absolute on the chromosome)
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(kwargs.pop('max_iter', 100))
		arg_list.append(kwargs.pop('ftol', 1e-6))
		arg_list.append(kwargs.pop('xtol', 1e-6))
		self._orig_init(*arg_list)
	nlopt_bobyqa._orig_init = nlopt_bobyqa.__init__
	nlopt_bobyqa.__init__ = _nlopt_bobyqa_ctor

	def _nlopt_sbplx_ctor(self,**kwargs):
		"""
		Constructs a Subplex (a variant of Nelder-Mead that uses Nelder-Mead on a sequence of subspaces) (NLOPT)
	
		USAGE: algorithm.nlopt_sbplx(max_iter = 100, ftol = 1e-6, xtol = 1e-6);
	
		* max_iter: stop-criteria (number of iterations)
		* ftol: stop-criteria (absolute on the obj-fun)
		* xtol: stop-criteria (absolute on the chromosome)
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(kwargs.pop('max_iter', 100))
		arg_list.append(kwargs.pop('ftol', 1e-6))
		arg_list.append(kwargs.pop('xtol', 1e-6))
		self._orig_init(*arg_list)
	nlopt_sbplx._orig_init = nlopt_sbplx.__init__
	nlopt_sbplx.__init__ = _nlopt_sbplx_ctor

	def _nlopt_cobyla_ctor(self,**kwargs):
		"""
		Constructs a Constrained Optimization BY Linear Approximation (COBYLA) algorithm (NLOPT)
	
		USAGE: algorithm.nlopt_cobyla(max_iter = 100, ftol = 1e-6, xtol = 1e-6)
	
		* max_iter: stop-criteria (number of iterations)
		* ftol: stop-criteria (absolute on the obj-fun)
		* xtol: stop-criteria (absolute on the chromosome)
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(kwargs.pop('max_iter', 100))
		arg_list.append(kwargs.pop('ftol', 1e-6))
		arg_list.append(kwargs.pop('xtol', 1e-6))
		self._orig_init(*arg_list)
	nlopt_cobyla._orig_init = nlopt_cobyla.__init__
	nlopt_cobyla.__init__ = _nlopt_cobyla_ctor

	def _nlopt_mma_ctor(self,**kwargs):
		"""
		Constructs a Method of Moving Asymptotes (MMA) algorithm (NLOPT)

		USAGE: algorithm.nlopt_mma(max_iter = 100, ftol = 1e-6, xtol = 1e-6)

		* max_iter: stop-criteria (number of iterations)
		* ftol: stop-criteria (absolute on the obj-fun)
		* xtol: stop-criteria (absolute on the chromosome)
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(kwargs.pop('max_iter', 100))
		arg_list.append(kwargs.pop('ftol', 1e-6))
		arg_list.append(kwargs.pop('xtol', 1e-6))
		self._orig_init(*arg_list)
	nlopt_mma._orig_init = nlopt_mma.__init__
	nlopt_mma.__init__ = _nlopt_mma_ctor

	def _nlopt_slsqp_ctor(self,**kwargs):
		"""
		Constructs a Sequential Least SQuares Programming algorithm (SLSQP) algorithm (NLOPT)

		USAGE: algorithm.nlopt_slsqp(max_iter = 100, ftol = 1e-6, xtol = 1e-6)

		* max_iter: stop-criteria (number of iterations)
		* ftol: stop-criteria (absolute on the obj-fun)
		* xtol: stop-criteria (absolute on the chromosome)
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(kwargs.pop('max_iter', 100))
		arg_list.append(kwargs.pop('ftol', 1e-6))
		arg_list.append(kwargs.pop('xtol', 1e-6))
		self._orig_init(*arg_list)
	nlopt_slsqp._orig_init = nlopt_slsqp.__init__
	nlopt_slsqp.__init__ = _nlopt_slsqp_ctor

#GSL algorithms (only if PyGMO has been compiled with gsl option activated)
if "gsl" in str(_get_algorithm_list()):
	def _gsl_bfgs_ctor(self,**kwargs):
		"""
		Constructs a BFGS Algorithm (GSL)
	
		USAGE: algorithm.gsl_bfgs(max_iter = 100, step_size = 1e-8, tol = 1e-8, grad_step_size = 0.01, grad_tol = 0.0001);
	
		* max_iter: maximum number of iterations
		* step_size: size of the first trial step.
		* tol: accuracy of the line minimisation.
		* grad_step_size: step size for the numerical computation of the gradient.
		* grad_tol: tolerance when testing the norm of the gradient as stopping criterion.
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(kwargs.pop('max_iter', 100))
		arg_list.append(kwargs.pop('grad_tol', 1e-4))
		arg_list.append(kwargs.pop('grad_step_size', 1e-2))
		arg_list.append(kwargs.pop('tol', 1e-8))
		arg_list.append(kwargs.pop('step_size', 1e-8))
		self._orig_init(*arg_list)
	gsl_bfgs._orig_init = gsl_bfgs.__init__
	gsl_bfgs.__init__ = _gsl_bfgs_ctor

	def _gsl_bfgs2_ctor(self,**kwargs):
		"""
		Constructs a BFGS2 Algorithm (GSL)

		NOTE: in GSL, BFGS2 is a more efficient version of BFGS
	
		USAGE: algorithm.gsl_bfgs2(max_iter = 100, step_size = 1e-8, tol = 1e-8, grad_step_size = 0.01, grad_tol = 0.0001);
	
		* max_iter: maximum number of iterations
		* step_size: size of the first trial step.
		* tol: accuracy of the line minimisation.
		* grad_step_size: step size for the numerical computation of the gradient.
		* grad_tol: tolerance when testing the norm of the gradient as stopping criterion.
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(kwargs.pop('max_iter', 100))
		arg_list.append(kwargs.pop('grad_tol', 1e-4))
		arg_list.append(kwargs.pop('grad_step_size', 1e-2))
		arg_list.append(kwargs.pop('tol', 1e-8))
		arg_list.append(kwargs.pop('step_size', 1e-8))
		self._orig_init(*arg_list)
	gsl_bfgs2._orig_init = gsl_bfgs2.__init__
	gsl_bfgs2.__init__ = _gsl_bfgs2_ctor

	def _gsl_fr_ctor(self,**kwargs):
		"""
		Constructs a Fletcher-Reeves conjugate gradient (GSL)
	
		USAGE: algorithm.gsl_fr(max_iter = 100, step_size = 1e-8, tol = 1e-8, grad_step_size = 0.01, grad_tol = 0.0001);
	
		* max_iter: maximum number of iterations
		* step_size: size of the first trial step.
		* tol: accuracy of the line minimisation.
		* grad_step_size: step size for the numerical computation of the gradient.
		* grad_tol: tolerance when testing the norm of the gradient as stopping criterion.
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(kwargs.pop('max_iter', 100))
		arg_list.append(kwargs.pop('grad_tol', 1e-4))
		arg_list.append(kwargs.pop('grad_step_size', 1e-2))
		arg_list.append(kwargs.pop('tol', 1e-8))
		arg_list.append(kwargs.pop('step_size', 1e-8))
		self._orig_init(*arg_list)
	gsl_fr._orig_init = gsl_fr.__init__
	gsl_fr.__init__ = _gsl_fr_ctor

	def _gsl_pr_ctor(self,**kwargs):
		"""
		Constructs a Polak-Ribiere conjugate gradient (GSL)

		USAGE: algorithm.gsl_pr2(max_iter = 100, step_size = 1e-8, tol = 1e-8, grad_step_size = 0.01, grad_tol = 0.0001);
	
		* max_iter: maximum number of iterations
		* step_size: size of the first trial step.
		* tol: accuracy of the line minimisation.
		* grad_step_size: step size for the numerical computation of the gradient.
		* grad_tol: tolerance when testing the norm of the gradient as stopping criterion.
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(kwargs.pop('max_iter', 100))
		arg_list.append(kwargs.pop('grad_tol', 1e-4))
		arg_list.append(kwargs.pop('grad_step_size', 1e-2))
		arg_list.append(kwargs.pop('tol', 1e-8))
		arg_list.append(kwargs.pop('step_size', 1e-8))
		self._orig_init(*arg_list)
	gsl_pr._orig_init = gsl_pr.__init__
	gsl_pr.__init__ = _gsl_pr_ctor

	def _gsl_nm_ctor(self,**kwargs):
		"""
		Constructs a Nelder-Mead Algorithm (GSL)

		USAGE: algorithm.gsl_nm(max_iter = 100, step_size = 1e-8, tol = 1e-8);
	
		* max_iter: maximum number of iterations
		* step_size: size of the first trial step.
		* tol: accuracy of the line minimisation.
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(kwargs.pop('max_iter', 100))
		arg_list.append(kwargs.pop('tol', 1e-8))
		arg_list.append(kwargs.pop('step_size', 1e-8))
		self._orig_init(*arg_list)
	gsl_nm._orig_init = gsl_nm.__init__
	gsl_nm.__init__ = _gsl_nm_ctor

	def _gsl_nm2_ctor(self,**kwargs):
		"""
		Constructs a Nelder-Mead algorithm (Variant2) (GSL)

		USAGE: algorithm.gsl_nm2(max_iter = 100, step_size = 1e-8, tol = 1e-8);
	
		* max_iter: maximum number of iterations
		* step_size: size of the first trial step.
		* tol: accuracy of the line minimisation.
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(kwargs.pop('max_iter', 100))
		arg_list.append(kwargs.pop('tol', 1e-8))
		arg_list.append(kwargs.pop('step_size', 1e-8))
		self._orig_init(*arg_list)
	gsl_nm2._orig_init = gsl_nm2.__init__
	gsl_nm2.__init__ = _gsl_nm2_ctor

	def _gsl_nm2rand_ctor(self,**kwargs):
		"""
		Constructs a Nelder-Mead algorithm (Variant2 + randomly oriented initial simplex) (GSL)

		USAGE: algorithm.gsl_nm2rand(max_iter = 100, step_size = 1e-8, tol = 1e-8);
	
		* max_iter: maximum number of iterations
		* step_size: size of the first trial step.
		* tol: accuracy of the line minimisation.
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(kwargs.pop('max_iter', 100))
		arg_list.append(kwargs.pop('tol', 1e-8))
		arg_list.append(kwargs.pop('step_size', 1e-8))
		self._orig_init(*arg_list)
	gsl_nm2rand._orig_init = gsl_nm2rand.__init__
	gsl_nm2rand.__init__ = _gsl_nm2rand_ctor

#IPOPT algorithm (only if PyGMO has been compiled with the ipopt option activated)
if "ipopt" in str(_get_algorithm_list()):
	def _ipopt_ctor(self,**kwargs):
		"""
		Constructs an Interior Point OPTimization Algorithm (IPOPT)
	
		USAGE: algorithm.ipopt(major_iter = 100, constr_viol_tol = 1e-08, dual_inf_tol = 1e-08, compl_inf_tol = 1e-08, screen_output = False);
	
		* major_iter: Maximum number of major iterations
		* constr_viol_tol: Constraint violation tolerance
		* dual_inf_tol: Dual infeasibility tolerance
		* compl_inf_tol: Complementary feasibility tolerance
		* screen_output: Activates output on screen
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(kwargs.pop('major_iter', 100))
		arg_list.append(kwargs.pop('constr_viol_tol', 1e-8))
		arg_list.append(kwargs.pop('dual_inf_tol', 1e-8))
		arg_list.append(kwargs.pop('compl_inf_tol', 1e-8))
		self._orig_init(*arg_list)
		self.screen_output = kwargs.pop('screen_output', False)
	ipopt._orig_init = ipopt.__init__
	ipopt.__init__ = _ipopt_ctor

#SNOPT algorithm (only if PyGMO has been compiled with the snopt option activated)
if "snopt" in str(_get_algorithm_list()):
	def _snopt_ctor(self,**kwargs):
		"""
		Constructs SNOPT Algorithm
	
		USAGE: algorithm.snopt(major_iter = 100, feas_tol = 1e-6, opt_tol = 1e-6, screen_output = False);
	
		* major_iter: Maximum number of major iterations
		* feas_tol: Feasibility tolerance
		* opt_tol: Optimality tolerance
		* screen_output: Activates output on screen
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(kwargs.pop('major_iter', 100))
		arg_list.append(kwargs.pop('feas_tol', 1e-8))
		arg_list.append(kwargs.pop('opt_tol', 1e-8))
		self._orig_init(*arg_list)
		self.screen_output = kwargs.pop('screen_output', False)
	snopt._orig_init = snopt.__init__
	snopt.__init__ = _snopt_ctor


	
	

