# -*- coding: utf-8 -*-
from _algorithm import *
from _algorithm import _base
from _base import base
from _example import py_example
from _cmaes import py_cmaes
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
	except ImportError as e:
		algorithm_list = [algorithm.__dict__[n] for n in filter(lambda n: not n.startswith('_') and not n == 'base' and not n.startswith('scipy'),dir(algorithm))]
	return algorithm_list
	


# Redefining the constructors of all algorithms to obtain good documentation and to allow kwargs
def _de_ctor(self, gen=100, f=0.8, cr=0.9, variant=2, ftol=1e-6, xtol=1e-6, screen_output = False):
	"""
	Constructs a Differential Evolution algorithm:
	
	USAGE: algorithm.de(gen=1, f=0.5, cr=0.9, variant=2, ftol=1e-6, xtol=1e-6, screen_output = False)
	
	* gen: number of generations
	* f: weighting factor in [0,1] (if -1 self-adptation is used)
	* cr: crossover in [0,1] (if -1 self-adptation is used)
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
	* ftol stop criteria on f
	* xtol stop criteria on x
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(gen)
	arg_list.append(f)
	arg_list.append(cr)
	arg_list.append(variant)
	arg_list.append(ftol)
	arg_list.append(xtol)
	self._orig_init(*arg_list)
	self.screen_output = screen_output
de._orig_init = de.__init__
de.__init__ = _de_ctor



def _jde_ctor(self, gen=100, variant=2, variant_adptv=1, ftol=1e-6, xtol=1e-6, memory=False, screen_output = False):
	"""
	Constructs a jDE algorithm (self-adaptive DE)
	
	REF: "Self-adaptive differential evolution algorithm in constrained real-parameter optimization"
	J Brest, V Zumer, MS Maucec - Evolutionary Computation, 2006. 
	http://dsp.szu.edu.cn/DSP2006/research/publication/yan/WebEdit/UploadFile/Self-adaptive%20Differential%20Evolution%20Algorithm%20for%20Constrained%20Real-Parameter%20Optimization.pdf
	
	USAGE: algorithm.jde(gen=100, variant=2, variant_adptv=1, ftol=1e-6, xtol=1e-6, memory = False, screen_output = False)
	
	* gen: number of generations
	* variant: algoritmic variant to use (one of [1 .. 18])
		1. best/1/exp				2. rand/1/exp
		3. rand-to-best/1/exp		4. best/2/exp
		5. rand/2/exp				6. best/1/bin
		7. rand/1/bin				8. rand-to-best/1/bin
		9. best/2/bin				10. rand/2/bin
		11. best/3/exp				12. best/3/bin
		13. rand/3/exp				14. rand/3/bin
		15. rand-to-current/2/exp		16. rand-to-current/2/bin
		17. rand-to-best-and-current/2/exp	18. rand-to-best-and-current/2/bin
	* variant_adptv: adaptive scheme to use (one of [1..2])
		1. random param mutation		2. param mutation follows rand/3 scheme
	* ftol: stop criteria on f
	* xtol: stop criteria on x
	* memory: if True the algorithm internal state is saved and used for the next call
	* screen_output: activates screen output of the algorithm (do not use in archipealgo, otherwise the screen will be flooded with 
	*				 different island outputs)
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(gen)
	arg_list.append(variant)
	arg_list.append(variant_adptv)
	arg_list.append(ftol)
	arg_list.append(xtol)
	arg_list.append(memory)	
	self._orig_init(*arg_list)
	self.screen_output = screen_output
jde._orig_init = jde.__init__
jde.__init__ = _jde_ctor

def _de_1220_ctor(self, gen=100, variant_adptv=1, allowed_variants = [1,2,3,4,5,6,7,8,9,10], memory = False, ftol=1e-6, xtol=1e-6, screen_output = False):
	"""
	Constructs a Differential Evolution algorithm (our own brew). Self adaptation on F, CR and mutation variant.:
	
	USAGE: algorithm.de_1220(gen=100, variant_adptv=1, allowed_variants = [i for i in range(1,19)], memory = False, ftol=1e-6, xtol=1e-6, screen_output = False)
	
	* gen: number of generations
	* variant_adptv: adaptiv scheme to use (one of [1..2])
		1. random param mutation		2. param mutation follows relative DE scheme
	* allowed_variants : a list of the algoritmic variants to mix and self-adapt. Allowed variants are ...
		1. best/1/exp				2. rand/1/exp
		3. rand-to-best/1/exp			4. best/2/exp
		5. rand/2/exp				6. best/1/bin
		7. rand/1/bin				8. rand-to-best/1/bin
		9. best/2/bin				10. rand/2/bin
		11. best/3/exp				12. best/3/bin
		13. rand/3/exp				14. rand/3/bin
		15. rand-to-current/2/exp		16. rand-to-current/2/bin
		17. rand-to-best-and-current/2/exp	18. rand-to-best-and-current/2/bin
	* ftol: stop criteria on f
	* xtol: stop criteria on x
	* memory: if True the algorithm internal state is saved and used for the next call
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(gen)
	arg_list.append(variant_adptv)
	arg_list.append(allowed_variants)
	arg_list.append(memory)
	arg_list.append(ftol)
	arg_list.append(xtol)
	self._orig_init(*arg_list)
	self.screen_output = screen_output
de_1220._orig_init = de_1220.__init__
de_1220.__init__ = _de_1220_ctor

def _mde_pbx_ctor(self, gen=100, qperc=0.15, nexp=1.5, ftol=1e-6, xtol=1e-6, screen_output = False):
	"""
	Constructs a mde_pbx algorithm (self-adaptive DE)
	
	REF: "An Adaptive Differential Evolution Algorithm With Novel Mutation and Crossover
	Strategies for Global Numerical Optimization" - IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS?PART B: CYBERNETICS, VOL. 42, NO. 2, APRIL 20 

	
	USAGE: algorithm.mde_pbx(gen=100, qperc=0.15, nexp=1.5, ftol=1e-6, xtol=1e-6, screen_output = False)
	
	* gen: number of generations
	* qperc: percentage of population to choose the best vector
	* nexp: exponent for the powermean
	* ftol: stop criteria on f
	* xtol: stop criteria on x
	* screen_output: activates screen output of the algorithm (do not use in archipealgo, otherwise the screen will be flooded with 
	* 		 different island outputs)
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(gen)
	arg_list.append(qperc)
	arg_list.append(nexp)
	arg_list.append(ftol)
	arg_list.append(xtol)	
	self._orig_init(*arg_list)
	self.screen_output = screen_output
mde_pbx._orig_init = mde_pbx.__init__
mde_pbx.__init__ = _mde_pbx_ctor


def _pso_ctor(self, gen=1, omega = 0.7298, eta1 = 2.05, eta2 = 2.05, vcoeff = 0.5, variant = 5, neighb_type = 2, neighb_param = 4):
	"""
	Constructs a Particle Swarm Optimization (steady-state). The position update is applied
	immediately after the velocity update
	
	REF (for variants 5-6): http://cswww.essex.ac.uk/staff/rpoli/papers/PoliKennedyBlackwellSI2007.pdf
	
	REF (for variants 1-4): Kennedy, J.; Eberhart, R. (1995). "Particle Swarm Optimization". Proceedings of IEEE International Conference on Neural Networks. IV. pp. 1942?1948.
	
	USAGE: algorithm.pso(gen=1, omega = 0.7298, eta1 = 2.05, eta2 = 2.05, vcoeff = 0.5, variant = 5, neighb_type = 2, neighb_param = 4)

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
	arg_list.append(gen)
	arg_list.append(omega)
	arg_list.append(eta1)
	arg_list.append(eta2)
	arg_list.append(vcoeff)
	arg_list.append(variant)
	arg_list.append(neighb_type)
	arg_list.append(neighb_param)	
	self._orig_init(*arg_list)
pso._orig_init = pso.__init__
pso.__init__ = _pso_ctor

def _pso_gen_ctor(self, gen=1, omega = 0.7298, eta1 = 2.05, eta2 = 2.05, vcoeff = 0.5, variant = 5, neighb_type = 2, neighb_param = 4, use_racing = False, max_fevals = -1):
	"""
	Constructs a Particle Swarm Optimization (generational). The position update is applied
	only at the end of an entire loop over the population (swarm). Use this version for stochastic problems.
	
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
	* use_racing: Whether to use racing
	* max_fevals: When specified other than -1, this serve as another termination condition -- maximium number of objective function evaluations
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(gen)
	arg_list.append(omega)
	arg_list.append(eta1)
	arg_list.append(eta2)
	arg_list.append(vcoeff)
	arg_list.append(variant)
	arg_list.append(neighb_type)
	arg_list.append(neighb_param)	
	arg_list.append(use_racing)
	if max_fevals > 0:
		arg_list.append(max_fevals)
	self._orig_init(*arg_list)
pso_gen._orig_init = pso_gen.__init__
pso_gen.__init__ = _pso_gen_ctor

def _sga_ctor(self, gen=1, cr=.95, m=.02, elitism=1, mutation=sga.mutation.GAUSSIAN, width = 0.1, selection=sga.selection.ROULETTE, crossover=sga.crossover.EXPONENTIAL):
	"""
	Constructs a Simple Genetic Algorithm (generational)
	
	USAGE: algorithm.sga(self, gen=1, cr=.95, m=.02, elitism=1, mutation=sga.mutation.GAUSSIAN, width = 0.1, selection=sga.selection.ROULETTE, crossover=sga.crossover.EXPONENTIAL)
  
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
	arg_list.append(gen)
	arg_list.append(cr)
	arg_list.append(m)
	arg_list.append(elitism)
	arg_list.append(mutation)
	arg_list.append(width)
	arg_list.append(selection)
	arg_list.append(crossover)	
	self._orig_init(*arg_list)
sga._orig_init = sga.__init__
sga.__init__ = _sga_ctor

def _vega_ctor(self, gen=1, cr=.95, m=.02, elitism=1, mutation=sga.mutation.GAUSSIAN, width = 0.1, selection=sga.selection.ROULETTE, crossover=sga.crossover.EXPONENTIAL):
	"""
	Constructs a Vector evaluated genetic algorithm
	
	USAGE: algorithm.vega(self, gen=1, cr=.95, m=.02, elitism=1, mutation=sga.mutation.GAUSSIAN, width = 0.1, crossover=sga.crossover.EXPONENTIAL)
  
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
	arg_list.append(gen)
	arg_list.append(cr)
	arg_list.append(m)
	arg_list.append(elitism)
	arg_list.append(mutation)
	arg_list.append(width)
	arg_list.append(crossover)	
	self._orig_init(*arg_list)
vega._orig_init = vega.__init__
vega.__init__ = _vega_ctor

def _nsga_II_ctor(self, gen=100, cr = 0.95, eta_c = 10, m = 0.01, eta_m = 10):
	"""
	Constructs a Non-dominated Sorting Genetic Algorithm (NSGA_II)
	
	USAGE: algorithm.nsga_II(self, gen=100, cr = 0.95, eta_c = 10, m = 0.01, eta_m = 10)
  
	* gen: number of generations
	* cr: crossover factor [0,1[
	* eta_c: Distribution index for crossover
	* m: mutation probability [0,1]
	* eta_m: Distribution index for mutation
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(gen)
	arg_list.append(cr)
	arg_list.append(eta_c)
	arg_list.append(m)
	arg_list.append(eta_m)
	self._orig_init(*arg_list)
nsga_II._orig_init = nsga_II.__init__
nsga_II.__init__ = _nsga_II_ctor


_algorithm.pade.RANDOM = _algorithm._weight_generation.RANDOM
_algorithm.pade.GRID = _algorithm._weight_generation.GRID
_algorithm.pade.LOW_DISCREPANCY = _algorithm._weight_generation.LOW_DISCREPANCY
from PyGMO.problem import decompose
def _pade_ctor(self, gen=10, max_parallelism = 1, decomposition = decompose.WEIGHTED, solver = jde(10), T = 8, weights = pade.RANDOM, z = None):
	"""
	Constructs a Parallel Decomposition Algorithm (PaDe).
	
	For each element of the population a different single objective problem is generated using a decomposition method.
	Those single-objective problems are thus solved in an island model.
	At the end of the evolution the population is set as the best individual in each single-objective island.
	This algorithm, original with PaGMO, builds upon the MOEA/D framework
	
	USAGE: algorithm.pade(self, gen=10, max_parallelism = 1, decomposition = decompose.WEIGHTED, solver = jde(10), T = 8, weights = pade.RANDOM, z = None)

	* gen: number of generations
	* max_parallelism: the maximum number of single-objective problems to solve at the same time
	* solver: the algorithm to use to solve the single-objective problems
	* T: the size of the population on each subproblem (must be an even number)
	* decomposition = the decomposition method to use (Weighted, Tchebycheff or BI)
	* weights: the weight generation method
	* z: the reference point (used with Tchebycheff and BI decomposition methods)
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(gen)
	arg_list.append(max_parallelism)
	arg_list.append(decomposition)
	arg_list.append(solver)
	arg_list.append(T)
	arg_list.append(weights)
	if z != None:
		arg_list.append(z)
	self._orig_init(*arg_list)
pade._orig_init = pade.__init__
pade.__init__ = _pade_ctor
del decompose

def _sa_corana_ctor(self, iter = 10000, Ts = 10, Tf = .1, steps = 1, bin_size = 20, range = 1):
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
	arg_list.append(iter)
	arg_list.append(Ts)
	arg_list.append(Tf)
	arg_list.append(steps)
	arg_list.append(bin_size)
	arg_list.append(range)
	self._orig_init(*arg_list)
sa_corana._orig_init = sa_corana.__init__
sa_corana.__init__ = _sa_corana_ctor

def _bee_colony_ctor(self, gen = 100, limit = 20):
	"""
	Constructs an Artificial Bee Colony Algorithm
	
	USAGE: algorithm.bee_colony(gen = 100, limit = 20)
	
	* gen: number of 'generations' (each generation 2*NP function evaluations
		are made where NP is the population size)
	* limit: number of tries after which a source of food is dropped if not improved
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(gen)
	arg_list.append(limit)
	self._orig_init(*arg_list)
bee_colony._orig_init = bee_colony.__init__
bee_colony.__init__ = _bee_colony_ctor

def _sea_ctor(self, gen = 100, limit = 20):
	"""
	Constructs a simple (N+1)-EA: A Simple Evolutionary Algorithm

	USAGE: algorithm.ea(gen = 1)
	SEE : Oliveto, Pietro S., Jun He, and Xin Yao.
	"Time complexity of evolutionary algorithms for combinatorial optimization: A decade of results."
	International Journal of Automation and Computing 4.3 (2007): 281-293.

	* gen: number of 'generations' (each generation is one function evaluation)
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(gen)
	self._orig_init(*arg_list)
sea._orig_init = sea.__init__
sea.__init__ = _sea_ctor

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

def _ms_ctor(self, algorithm = _algorithm.de(), iter = 1):
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
	arg_list.append(algorithm)
	arg_list.append(iter)
	self._orig_init(*arg_list)
ms._orig_init = ms.__init__
ms.__init__ = _ms_ctor

def _cs_ctor(self, max_eval = 1, stop_range = 0.01, start_range = 0.1, reduction_coeff = 0.5):
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
	arg_list.append(max_eval)
	arg_list.append(stop_range)
	arg_list.append(start_range)
	arg_list.append(reduction_coeff)
	self._orig_init(*arg_list)
cs._orig_init = cs.__init__
cs.__init__ = _cs_ctor

def _mbh_ctor(self, algorithm = _algorithm.cs(), stop = 5, perturb = 5e-2, screen_output = False):
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
	* screen_output: activates screen output of the algorithm (do not use in archipealgo, otherwise the screen will be flooded with 
	* 		 different island outputs)
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(algorithm)
	arg_list.append(stop)
	arg_list.append(perturb)
	self._orig_init(*arg_list)
	self.screen_output = screen_output
mbh._orig_init = mbh.__init__
mbh.__init__ = _mbh_ctor

def _ihs_ctor(self, iter = 100, hmcr = 0.85, par_min = 0.35, par_max = 0.99, bw_min = 1E-5, bw_max = 1):
	"""
	Constructs an Improved Harmony Search Algorithm
	
	USAGE: algorithm.ihs(iter = 100, hmcr = 0.85, par_min = 0.35, par_max = 0.99, bw_min = 1E-5, bw_max = 1);
	
	* iter: number of iterations (improvisations)
	* hmcr: rate of choosing from memory (in ]0,1[)
	* par_min: minimum pitch adjustment rate (in ]0,1[)
	* par_max: maximum pitch adjustment rate (in ]0,1[, > par_min) 
	* bw_min: minimum distance bandwidth 
	* bw_max: maximum distance bandwidth (> bw_min)
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(iter)
	arg_list.append(hmcr)
	arg_list.append(par_min)
	arg_list.append(par_max)
	arg_list.append(bw_min)
	arg_list.append(bw_max)
	self._orig_init(*arg_list)
ihs._orig_init = ihs.__init__
ihs.__init__ = _ihs_ctor

def _cmaes_ctor(self, gen = 500, cc = -1, cs = -1, c1 = -1, cmu = -1, sigma0=0.5, ftol = 1e-6, xtol = 1e-6, memory = False, screen_output = False):
	"""
	Constructs a Covariance Matrix Adaptation Evolutionary Strategy (C++)

	USAGE: algorithm.cmaes(gen = 500, cc = -1, cs = -1, c1 = -1, cmu = -1, sigma0=0.5, ftol = 1e-6, xtol = 1e-6, memory = False, screen_output = False)

	NOTE: In our variant of the algorithm, particle memory is used to extract the elite and reinsertion
	is made aggressively ..... getting rid of the worst guy). Also, the bounds of the problem
	are enforced, as to allow PaGMO machinery to work. Fine control on each iteration can be achieved
	by calling the algo with memory=True and gen=1

	* gen: number of generations
	* cc: time constant for C cumulation (in [0,1]) if -1 automatic values are set
	* cs: time constant for sigma cumulation (in [0,1]) if -1 automatic values are set
	* c1: learning rate for rank-1 update (in [0,1]) if -1 automatic values are set
	* cmu: learning rate for rank-mu update (in [0,1]) if -1 automatic values are set
	* sigma0: starting step (std)
	* xtol: stopping criteria on the x tolerance
	* ftol: stopping criteria on the f tolerance
	* memory: if True the algorithm internal state is saved and used for the next call	
	* screen_output: activates screen output of the algorithm (do not use in archipealgo, otherwise the screen will be flooded with 
	* 		 different island outputs)
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(gen)
	arg_list.append(cc)
	arg_list.append(cs)
	arg_list.append(c1)
	arg_list.append(cmu)
	arg_list.append(sigma0)
	arg_list.append(ftol)
	arg_list.append(xtol)
	arg_list.append(memory)
	self._orig_init(*arg_list)
	self.screen_output = screen_output
cmaes._orig_init = cmaes.__init__
cmaes.__init__ = _cmaes_ctor

def _monte_carlo_ctor(self, iter = 10000):
	"""
	Constructs a Monte Carlo Algorithm
	
	USAGE: algorithm.monte_carlo(iter = 10000)
	
	NOTE: At the end of each iteration, the randomly generated 
		point substitutes the worst in the population if better
	
	* iter: number of Monte Carlo runs
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(iter)
	self._orig_init(*arg_list)
monte_carlo._orig_init = monte_carlo.__init__
monte_carlo.__init__ = _monte_carlo_ctor

#NLOPT algorithms (only if PyGMO has been compiled woth nlopt option activated)
if "nlopt" in str(_get_algorithm_list()):
	def _nlopt_bobyqa_ctor(self, max_iter = 100, ftol = 1e-6, xtol = 1e-6):
		"""
		Constructs a BOBYQA algorithm (Bound Optimization BY Quadratic Approximation) (NLOPT)
	
		USAGE: algorithm.nlopt_bobyqa(max_iter = 100, ftol = 1e-6, xtol = 1e-6)
	
		* max_iter: stop-criteria (number of iterations)
		* ftol: stop-criteria (absolute on the obj-fun)
		* xtol: stop-criteria (absolute on the chromosome)
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(max_iter)
		arg_list.append(ftol)
		arg_list.append(xtol)
		self._orig_init(*arg_list)
	nlopt_bobyqa._orig_init = nlopt_bobyqa.__init__
	nlopt_bobyqa.__init__ = _nlopt_bobyqa_ctor

	def _nlopt_sbplx_ctor(self, max_iter = 100, ftol = 1e-6, xtol = 1e-6):
		"""
		Constructs a Subplex (a variant of Nelder-Mead that uses Nelder-Mead on a sequence of subspaces) (NLOPT)
	
		USAGE: algorithm.nlopt_sbplx(max_iter = 100, ftol = 1e-6, xtol = 1e-6)
	
		* max_iter: stop-criteria (number of iterations)
		* ftol: stop-criteria (absolute on the obj-fun)
		* xtol: stop-criteria (absolute on the chromosome)
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(max_iter)
		arg_list.append(ftol)
		arg_list.append(xtol)
		self._orig_init(*arg_list)
	nlopt_sbplx._orig_init = nlopt_sbplx.__init__
	nlopt_sbplx.__init__ = _nlopt_sbplx_ctor

	def _nlopt_cobyla_ctor(self, max_iter = 100, ftol = 1e-6, xtol = 1e-6):
		"""
		Constructs a Constrained Optimization BY Linear Approximation (COBYLA) algorithm (NLOPT)
	
		USAGE: algorithm.nlopt_cobyla(max_iter = 100, ftol = 1e-6, xtol = 1e-6)
	
		* max_iter: stop-criteria (number of iterations)
		* ftol: stop-criteria (absolute on the obj-fun)
		* xtol: stop-criteria (absolute on the chromosome)
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(max_iter)
		arg_list.append(ftol)
		arg_list.append(xtol)
		self._orig_init(*arg_list)
	nlopt_cobyla._orig_init = nlopt_cobyla.__init__
	nlopt_cobyla.__init__ = _nlopt_cobyla_ctor

	def _nlopt_mma_ctor(self, max_iter = 100, ftol = 1e-6, xtol = 1e-6):
		"""
		Constructs a Method of Moving Asymptotes (MMA) algorithm (NLOPT)

		USAGE: algorithm.nlopt_mma(max_iter = 100, ftol = 1e-6, xtol = 1e-6)

		* max_iter: stop-criteria (number of iterations)
		* ftol: stop-criteria (absolute on the obj-fun)
		* xtol: stop-criteria (absolute on the chromosome)
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(max_iter)
		arg_list.append(ftol)
		arg_list.append(xtol)
		self._orig_init(*arg_list)
	nlopt_mma._orig_init = nlopt_mma.__init__
	nlopt_mma.__init__ = _nlopt_mma_ctor

	def _nlopt_auglag_ctor(self, aux_algo_id = 1, max_iter = 100, ftol = 1e-6, xtol = 1e-6, aux_max_iter = 100, aux_ftol = 1e-6, aux_xtol = 1e-6):
		"""
		Constructs an Augmented agrangian Algotihm (NLOPT)

		USAGE: algorithm.nlopt_mma(aux_algo_id = 1, max_iter = 100, ftol = 1e-6, xtol = 1e-6, aux_max_iter = 100, aux_ftol = 1e-6, aux_xtol = 1e-6)

		* aux_algo_id: auxiliary  optimizer id
			1: SBPLX
			2: COBYLA
			3: BOBYQA
			4: Low Storage BFGS
		* max_iter: stop-criteria (number of iterations)
		* ftol: stop-criteria (absolute on the obj-fun)
		* xtol: stop-criteria (absolute on the chromosome)
		* aux_max_iter: stop-criteria for the auxiliary optimizer (number of iterations)
		* aux_ftol: stop-criteria for the auxiliary optimizer (absolute on the obj-fun)
		* aux_xtol: stop-criteria for the auxiliary optimizer (absolute on the chromosome)
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(aux_algo_id)
		arg_list.append(max_iter)
		arg_list.append(ftol)
		arg_list.append(xtol)
		arg_list.append(aux_max_iter)
		arg_list.append(aux_ftol)
		arg_list.append(aux_xtol)
		self._orig_init(*arg_list)
	nlopt_auglag._orig_init = nlopt_auglag.__init__
	nlopt_auglag.__init__ = _nlopt_auglag_ctor

	def _nlopt_auglag_eq_ctor(self, aux_algo_id = 1, max_iter = 100, ftol = 1e-6, xtol = 1e-6, aux_max_iter = 100, aux_ftol = 1e-6, aux_xtol = 1e-6):
		"""
		Constructs an Augmented agrangian Algotihm (using penalties only for the equalities) (NLOPT)

		USAGE: algorithm.nlopt_auglag_eq(aux_algo_id = 1, max_iter = 100, ftol = 1e-6, xtol = 1e-6, aux_max_iter = 100, aux_ftol = 1e-6, aux_xtol = 1e-6)

		* aux_algo_id: auxiliary (local) optimizer id
			1: COBYLA
			2: MMA
		* max_iter: stop-criteria (number of iterations)
		* ftol: stop-criteria (absolute on the obj-fun)
		* xtol: stop-criteria (absolute on the chromosome)
		* aux_max_iter: stop-criteria for the auxiliary optimizer (number of iterations)
		* aux_ftol: stop-criteria for the auxiliary optimizer (absolute on the obj-fun)
		* aux_xtol: stop-criteria for the auxiliary optimizer (absolute on the chromosome)
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(aux_algo_id)
		arg_list.append(max_iter)
		arg_list.append(ftol)
		arg_list.append(xtol)
		arg_list.append(aux_max_iter)
		arg_list.append(aux_ftol)
		arg_list.append(aux_xtol)
		self._orig_init(*arg_list)
	nlopt_auglag_eq._orig_init = nlopt_auglag_eq.__init__
	nlopt_auglag_eq.__init__ = _nlopt_auglag_eq_ctor

	def _nlopt_slsqp_ctor(self, max_iter = 100, ftol = 1e-6, xtol = 1e-6):
		"""
		Constructs a Sequential Least SQuares Programming algorithm (SLSQP) algorithm (NLOPT)

		USAGE: algorithm.nlopt_slsqp(max_iter = 100, ftol = 1e-6, xtol = 1e-6)

		* max_iter: stop-criteria (number of iterations)
		* ftol: stop-criteria (absolute on the obj-fun)
		* xtol: stop-criteria (absolute on the chromosome)
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(max_iter)
		arg_list.append(ftol)
		arg_list.append(xtol)
		self._orig_init(*arg_list)
	nlopt_slsqp._orig_init = nlopt_slsqp.__init__
	nlopt_slsqp.__init__ = _nlopt_slsqp_ctor

#GSL algorithms (only if PyGMO has been compiled with gsl option activated)
if "gsl" in str(_get_algorithm_list()):
	def _gsl_bfgs_ctor(self, max_iter = 100, step_size = 1e-8, tol = 1e-8, grad_step_size = 0.01, grad_tol = 0.0001):
		"""
		Constructs a BFGS Algorithm (GSL)
	
		USAGE: algorithm.gsl_bfgs(max_iter = 100, step_size = 1e-8, tol = 1e-8, grad_step_size = 0.01, grad_tol = 0.0001)
	
		* max_iter: maximum number of iterations
		* step_size: size of the first trial step.
		* tol: accuracy of the line minimisation.
		* grad_step_size: step size for the numerical computation of the gradient.
		* grad_tol: tolerance when testing the norm of the gradient as stopping criterion.
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(max_iter)
		arg_list.append(grad_tol)
		arg_list.append(grad_step_size)
		arg_list.append(tol)
		arg_list.append(step_size)
		self._orig_init(*arg_list)
	gsl_bfgs._orig_init = gsl_bfgs.__init__
	gsl_bfgs.__init__ = _gsl_bfgs_ctor

	def _gsl_bfgs2_ctor(self, max_iter = 100, step_size = 1e-8, tol = 1e-8, grad_step_size = 0.01, grad_tol = 0.0001):
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
		arg_list.append(max_iter)
		arg_list.append(grad_tol)
		arg_list.append(grad_step_size)
		arg_list.append(tol)
		arg_list.append(step_size)
		self._orig_init(*arg_list)
	gsl_bfgs2._orig_init = gsl_bfgs2.__init__
	gsl_bfgs2.__init__ = _gsl_bfgs2_ctor

	def _gsl_fr_ctor(self, max_iter = 100, step_size = 1e-8, tol = 1e-8, grad_step_size = 0.01, grad_tol = 0.0001):
		"""
		Constructs a Fletcher-Reeves conjugate gradient (GSL)
	
		USAGE: algorithm.gsl_fr(max_iter = 100, step_size = 1e-8, tol = 1e-8, grad_step_size = 0.01, grad_tol = 0.0001)
	
		* max_iter: maximum number of iterations
		* step_size: size of the first trial step.
		* tol: accuracy of the line minimisation.
		* grad_step_size: step size for the numerical computation of the gradient.
		* grad_tol: tolerance when testing the norm of the gradient as stopping criterion.
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(max_iter)
		arg_list.append(grad_tol)
		arg_list.append(grad_step_size)
		arg_list.append(tol)
		arg_list.append(step_size)
		self._orig_init(*arg_list)
	gsl_fr._orig_init = gsl_fr.__init__
	gsl_fr.__init__ = _gsl_fr_ctor

	def _gsl_pr_ctor(self, max_iter = 100, step_size = 1e-8, tol = 1e-8, grad_step_size = 0.01, grad_tol = 0.0001):
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
		arg_list.append(max_iter)
		arg_list.append(grad_tol)
		arg_list.append(grad_step_size)
		arg_list.append(tol)
		arg_list.append(step_size)
		self._orig_init(*arg_list)
	gsl_pr._orig_init = gsl_pr.__init__
	gsl_pr.__init__ = _gsl_pr_ctor

	def _gsl_nm_ctor(self, max_iter = 100, step_size = 1e-8, tol = 1e-8):
		"""
		Constructs a Nelder-Mead Algorithm (GSL)

		USAGE: algorithm.gsl_nm(max_iter = 100, step_size = 1e-8, tol = 1e-8);
	
		* max_iter: maximum number of iterations
		* step_size: size of the first trial step.
		* tol: accuracy of the line minimisation.
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(max_iter)
		arg_list.append(tol)
		arg_list.append(step_size)
		self._orig_init(*arg_list)
	gsl_nm._orig_init = gsl_nm.__init__
	gsl_nm.__init__ = _gsl_nm_ctor

	def _gsl_nm2_ctor(self, max_iter = 100, step_size = 1e-8, tol = 1e-8):
		"""
		Constructs a Nelder-Mead algorithm (Variant2) (GSL)

		USAGE: algorithm.gsl_nm2(max_iter = 100, step_size = 1e-8, tol = 1e-8)
	
		* max_iter: maximum number of iterations
		* step_size: size of the first trial step.
		* tol: accuracy of the line minimisation.
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(max_iter)
		arg_list.append(tol)
		arg_list.append(step_size)
		self._orig_init(*arg_list)
	gsl_nm2._orig_init = gsl_nm2.__init__
	gsl_nm2.__init__ = _gsl_nm2_ctor

	def _gsl_nm2rand_ctor(self, max_iter = 100, step_size = 1e-8, tol = 1e-8):
		"""
		Constructs a Nelder-Mead algorithm (Variant2 + randomly oriented initial simplex) (GSL)

		USAGE: algorithm.gsl_nm2rand(max_iter = 100, step_size = 1e-8, tol = 1e-8);
	
		* max_iter: maximum number of iterations
		* step_size: size of the first trial step.
		* tol: accuracy of the line minimisation.
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(max_iter)
		arg_list.append(tol)
		arg_list.append(step_size)
		self._orig_init(*arg_list)
	gsl_nm2rand._orig_init = gsl_nm2rand.__init__
	gsl_nm2rand.__init__ = _gsl_nm2rand_ctor

#IPOPT algorithm (only if PyGMO has been compiled with the ipopt option activated)
if "ipopt" in str(_get_algorithm_list()):
	def _ipopt_ctor(self, max_iter = 100, constr_viol_tol = 1e-08, dual_inf_tol = 1e-08, compl_inf_tol = 1e-08, 
	nlp_scaling_method = True, obj_scaling_factor = 1.0, mu_init = 0.1, screen_output = False):
		"""
		Constructs an Interior Point OPTimization Algorithm (IPOPT)
	
		USAGE: algorithm.ipopt(major_iter = 100, constr_viol_tol = 1e-08, dual_inf_tol = 1e-08, compl_inf_tol = 1e-08, screen_output = False);
	
		* max_iter: Maximum number of major iterations
		* constr_viol_tol: Constraint violation tolerance
		* dual_inf_tol: Dual infeasibility tolerance
		* compl_inf_tol: Complementary feasibility tolerance
		* nlp_scaling_method Select if the "gradient-based" scaling of the  NLP should be used
		* obj_scaling_factor Scaling factor for the objective function.
		* mu_init Initial value for the barrier parameter.
		* screen_output: Activates output on screen
		"""
		# We set the defaults or the kwargs
		arg_list=[]
		arg_list.append(max_iter)
		arg_list.append(constr_viol_tol)
		arg_list.append(dual_inf_tol)
		arg_list.append(compl_inf_tol)
		arg_list.append(nlp_scaling_method)
		arg_list.append(obj_scaling_factor)
		arg_list.append(mu_init)
		self._orig_init(*arg_list)
		self.screen_output = screen_output
	ipopt._orig_init = ipopt.__init__
	ipopt.__init__ = _ipopt_ctor

#SNOPT algorithm (only if PyGMO has been compiled with the snopt option activated)
if "snopt" in str(_get_algorithm_list()):
	def _snopt_ctor(self,major_iter = 100, feas_tol = 1e-6, opt_tol = 1e-6, screen_output = False):
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
		arg_list.append(major_iter)
		arg_list.append(feas_tol)
		arg_list.append(opt_tol)
		self._orig_init(*arg_list)
		self.screen_output = screen_output
	snopt._orig_init = snopt.__init__
	snopt.__init__ = _snopt_ctor

