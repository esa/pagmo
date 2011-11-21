# -*- coding: iso-8859-1 -*-
from _base import base
from _base_stochastic import base_stochastic
from _problem import *
from _problem import _base
from _problem import _base_stochastic
from _example import py_example
from _example_stochastic import py_example_stochastic
from _pl2pl import py_pl2pl

# If GTOP database support is active import interplanetary trajectory problems

try:
	from _gtop import *
except:
	pass

# If GSL support is active import mit_sphere
try:
	from _mit_spheres import visualize as _visualize
	mit_spheres.visualize = _visualize
	def _mit_spheres_ctor(self, sample_size = 10, n_hidden = 10, ode_prec = 1E-3, seed = 0, symmetric = False, simulation_time = 50.0):
		"""
		Construct a Neurocontroller Evolution problem that seeks to drive three point masses to form a triangle
		This problem was used to design a contorller for the MIT SPHERES test bed on boear the ISS
		
		USAGE: problem.mit_spheres(sample_size = 10, n_hidden = 10, ode_prec = 1E-3, seed = 0, symmetric = False, simulation_time = 50.0):
		
		* sample_size: number of initial conditions the neurocontroller is tested from
		* n_hidden: number of hidden  for the feed-forward neural network
		* ode_prec: relative numerical precision of neurons the ODE integrator
		* seed: integer used as starting random seed to build the pseudorandom sequences used to generate the sample
		* symmetric: when True activates a Neural Network having symmetric weights (i.e. purely homogeneuos agents)
		* simulation_time: when True activates a Neural Network having symmetric weights (i.e. purely homogeneuos agents)
	"""
		
		# We construct the arg list for the original constructor exposed by boost_python
		arg_list=[]
		arg_list.append(sample_size)
		arg_list.append(n_hidden)
		arg_list.append(ode_prec)
		arg_list.append(seed)
		arg_list.append(symmetric)
		arg_list.append(simulation_time)
		self._orig_init(*arg_list)
	mit_spheres._orig_init = mit_spheres.__init__
	mit_spheres.__init__ = _mit_spheres_ctor
	
	from PyGMO import __version__
	__version__ = __version__ + "GTOP " + "GSL "

except:
	pass

	
#Creating the list of problems
def _get_problem_list():
	from PyGMO import problem
	return [problem.__dict__[n] for n in filter(lambda n: not n.startswith('_') and not n == 'base' and (issubclass(problem.__dict__[n],problem._base) or issubclass(problem.__dict__[n],problem._base_stochastic)),dir(problem))]

# Redefining the constructors of all problems to obtain good documentation and allowing kwargs
def _rastrigin_ctor(self,dim = 10):
	"""
	Constructs a Rastrigin problem (Box-Constrained Continuous Single-Objective)
	
	USAGE: problem.rastrigin(dim=10)
	
	* dim: problem dimension
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(dim)
	self._orig_init(*arg_list)
rastrigin._orig_init = rastrigin.__init__
rastrigin.__init__ = _rastrigin_ctor

def _rosenbrock_ctor(self,dim = 10):
	"""
	Constructs a Rosenbrock problem (Box-Constrained Continuous Single-Objective)
	
	USAGE: problem.rosenbrock(dim=10)
	
	* dim: problem dimension
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(dim)
	self._orig_init(*arg_list)
rosenbrock._orig_init = rosenbrock.__init__
rosenbrock.__init__ = _rosenbrock_ctor

def _ackley_ctor(self,dim = 10):
	"""
	Constructs a Ackley problem (Box-Constrained Continuous Single-Objective)
	
	USAGE: problem.ackley(dim=10)
	
	* dim: problem dimension
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(dim)
	self._orig_init(*arg_list)
ackley._orig_init = ackley.__init__
ackley.__init__ = _ackley_ctor

def _schwefel_ctor(self,dim = 10):
	"""
	Constructs a Schwefel problem (Box-Constrained Continuous Single-Objective)
	
	USAGE: problem.schwefel(dim=10)
	
	* dim: problem dimension
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(dim)
	self._orig_init(*arg_list)
schwefel._orig_init = schwefel.__init__
schwefel.__init__ = _schwefel_ctor

def _dejong_ctor(self,dim = 10):
	"""
	Constructs a De Jong problem (Box-Constrained Continuous Single-Objective)
	
	USAGE: problem.dejong(dim=10)
	
	* dim: problem dimension
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(dim)
	self._orig_init(*arg_list)
dejong._orig_init = dejong.__init__
dejong.__init__ = _dejong_ctor

def _griewank_ctor(self,dim = 10):
	"""
	Constructs a Griewank problem (Box-Constrained Continuous Single-Objective)
	
	USAGE: problem.griewank(dim=10)
	
	* dim: problem dimension
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(dim)
	self._orig_init(*arg_list)
griewank._orig_init = griewank.__init__
griewank.__init__ = _dejong_ctor

def _lennard_jones_ctor(self,n_atoms = 4):
	"""
	Constructs a Lennard-Jones problem (Box-Constrained Continuous Single-Objective)
	
	USAGE: problem.lennard_jones(n_atoms=4)
	
	* n_atoms: number of atoms
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(n_atoms)
	self._orig_init(*arg_list)
lennard_jones._orig_init = lennard_jones.__init__
lennard_jones.__init__ = _lennard_jones_ctor

def _branin_ctor(self):
	"""
	Constructs a Branin problem (Box-Constrained Continuous Single-Objective)
	
	USAGE: problem.branin()
	
	"""
	
branin._orig_init = branin.__init__
branin.__init__ = _branin_ctor

def _himmelblau_ctor(self):
	"""
	Constructs a Himmelblau problem (Box-Constrained Continuous Single-Objective)
	
	USAGE: problem.himmelblau()
	
	"""

himmelblau._orig_init = himmelblau.__init__
himmelblau.__init__ = _himmelblau_ctor

def _michalewicz_ctor(self,dim = 10):
	"""
	Constructs a Michalewicz problem (Box-Constrained Continuous Single-Objective)
	
	USAGE: problem.michalewicz(dim=5)
	
	NOTE: Minimum is -4.687 for dim=5 and -9.66 for dim = 10
	
	* dim: problem dimension
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(dim)
	self._orig_init(*arg_list)
michalewicz._orig_init = michalewicz.__init__
michalewicz.__init__ = _michalewicz_ctor

def _kur_ctor(self):
	"""
	Constructs a Kursawe's study problem (Box-Constrained Continuous Multi-Objective)
	
	USAGE: problem.kur()
	"""
	
kur._orig_init = kur.__init__
kur.__init__ = _kur_ctor

def _fon_ctor(self):
	"""
	Constructs a Fonseca and Fleming's study problem (Box-Constrained Continuous Multi-Objective)
	
	USAGE: problem.fon()
	"""
	
fon._orig_init = fon.__init__
fon.__init__ = _fon_ctor

def _pol_ctor(self):
	"""
	Constructs a Poloni's study study problem (Box-Constrained Continuous Multi-Objective)
	
	USAGE: problem.pol()
	"""
	
pol._orig_init = pol.__init__
pol.__init__ = _pol_ctor

def _sch_ctor(self):
	"""
	Constructs a Schaffer's study problem (Box-Constrained Continuous Multi-Objective)
	
	USAGE: problem.sch()
	"""
	
sch._orig_init = sch.__init__
sch.__init__ = _sch_ctor

def _zdt1_ctor(self):
	"""
	Constructs a ZDT1 problem (Box-Constrained Continuous Multi-Objective)
	
	USAGE: problem.zdt1()
	"""
	
zdt1._orig_init = zdt1.__init__
zdt1.__init__ = _zdt1_ctor

def _zdt2_ctor(self):
	"""
	Constructs a ZDT2 problem (Box-Constrained Continuous Multi-Objective)
	
	USAGE: problem.zdt2()
	"""
	
zdt2._orig_init = zdt2.__init__
zdt2.__init__ = _zdt2_ctor

def _zdt4_ctor(self):
	"""
	Constructs a ZDT4 problem (Box-Constrained Continuous Multi-Objective)
	
	USAGE: problem.zdt4()
	"""
	
zdt4._orig_init = zdt4.__init__
zdt4.__init__ = _zdt4_ctor

def _zdt6_ctor(self):
	"""
	Constructs a ZDT6 problem (Box-Constrained Continuous Multi-Objective)
	
	USAGE: problem.zdt6()
	"""
	
zdt6._orig_init = zdt6.__init__
zdt6.__init__ = _zdt6_ctor

def _luksan_vlcek_1_ctor(self,dim = 3):
	"""
	Constructs the first Luksan Vlcek problem (Constrained Continuous Single-Objective)
	
	NOTE: L. Luksan and J. Vlcek, "Sparse and Parially Separable Test Problems for Unconstrained and Equality Constrained Optimization"
	
	USAGE: problem.luksan_vlcek_1(dim=3)
	
	* dim: problem dimension
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(dim)
	self._orig_init(*arg_list)
luksan_vlcek_1._orig_init = luksan_vlcek_1.__init__
luksan_vlcek_1.__init__ = _luksan_vlcek_1_ctor

def _luksan_vlcek_2_ctor(self,dim = 16):
	"""
	Constructs the second Luksan Vlcek problem (Constrained Continuous Single-Objective)
	
	NOTE: L. Luksan and J. Vlcek, "Sparse and Parially Separable Test Problems for Unconstrained and Equality Constrained Optimization"
	
	USAGE: problem.luksan_vlcek_2(dim=16)
	
	* dim: problem dimension
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(dim)
	self._orig_init(*arg_list)
luksan_vlcek_2._orig_init = luksan_vlcek_2.__init__
luksan_vlcek_2.__init__ = _luksan_vlcek_2_ctor

def _luksan_vlcek_3_ctor(self,dim = 8):
	"""
	Constructs the third Luksan Vlcek problem (Constrained Continuous Single-Objective)
	
	NOTE: L. Luksan and J. Vlcek, "Sparse and Parially Separable Test Problems for Unconstrained and Equality Constrained Optimization"
	
	USAGE: problem.luksan_vlcek_3(dim=8)
	
	* dim: problem dimension
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(dim)
	self._orig_init(*arg_list)
luksan_vlcek_3._orig_init = luksan_vlcek_3.__init__
luksan_vlcek_3.__init__ = _luksan_vlcek_3_ctor

def _snopt_toyprob_ctor(self):
	"""
	Constructs SNOPT toy-problem (Box-Constrained Continuous Multi-Objective)
	
	USAGE: problem.snopt_toyprob()
	"""
	
snopt_toyprob._orig_init = snopt_toyprob.__init__
snopt_toyprob.__init__ = _snopt_toyprob_ctor

def _string_match_ctor(self,string = "Can we use it for space?"):
	"""
	Constructs a string-match problem (Box-Constrained Integer Single-Objective)
	
	NOTE: This is the problem of matching a string. Transcribed as an optimization problem
	
	USAGE: problem.string_match(string = "mah")
	
	* string: string to match
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(string)
	self._orig_init(*arg_list)
string_match._orig_init = string_match.__init__
string_match.__init__ = _string_match_ctor

def _golomb_ruler_ctor(self,order = 5, length=10):
	"""
	Constructs a Golomb Ruler problem (Constrained Integer Single-Objective)
	
	NOTE: see http://en.wikipedia.org/wiki/Golomb_ruler
	
	USAGE: problem.golomb_ruler(order = 5, length=10)
	
	* order: order of the Golomb ruler
	* length: length of the Golomb ruler
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(order)
	arg_list.append(length)
	self._orig_init(*arg_list)
golomb_ruler._orig_init = golomb_ruler.__init__
golomb_ruler.__init__ = _golomb_ruler_ctor

def _tsp_ctor(self,matrix = [[0,1,2],[1,0,5],[2,5,0]]):
	"""
	Constructs a Travelling Salesman problem (Constrained Integer Single-Objective)
	
	USAGE: problem.tsp(matrix = [0,1,2],[1,0,5],[2,5,0])
	
	* matrix: inter-city distances (symmetric matrix)
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(matrix)
	self._orig_init(*arg_list)
tsp._orig_init = tsp.__init__
tsp.__init__ = _tsp_ctor

def _knapsack_ctor(self,values = [1,2,3,4,5], weights = [10, 40, 30, 50, 20], max_weight = 100):
	"""
	Constructs a 0-1 Knapsack Problem (Constrained Integer Single-Objective)
	
	USAGE: problem.knapsack(values = [1,2,3,4,5], weights = [10, 40, 30, 50, 20], max_weight = 100)
	
	* values: raw array of values
	* weights: raw array of weights
	* max_weight: maximum weight
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(values)
	arg_list.append(weights)
	arg_list.append(max_weight)
	self._orig_init(*arg_list)
knapsack._orig_init = knapsack.__init__
knapsack.__init__ = _knapsack_ctor

def _inventory_ctor(self, weeks = 4, sample_size = 10, seed = 0):
	"""
	Constructs an Inventory Problem (Stochastic Objective Function)
	
	NOTE: see www2.isye.gatech.edu/people/faculty/Alex_Shapiro/SPbook.pdf
	
	USAGE: problem.inventory(weeks = 4, sample_size = 10, seed = 0):
	
	* week: dimension of the problem corresponding to the numer of weeks
	                 to plan the inventory for.
	* sample_size: dimension of the sample used to approximate the expected value
	* seed: integer used as starting random seed to build the 
	                 pseudorandom sequences used to generate the sample
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(weeks)
	arg_list.append(sample_size)
	arg_list.append(seed)
	self._orig_init(*arg_list)
inventory._orig_init = inventory.__init__
inventory.__init__ = _inventory_ctor


