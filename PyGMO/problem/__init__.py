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
except ImportError:
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

def _mo3d_plot(self, pop, a=40, comp=[0,1,2]):
	"""
	Generic plot-method for multi-objective optimization problems with more then 2 objectives

	USAGE: prob.plot(pop, comp[0,2,3])
	* pop: population of solutions to the problem
	* a: angle of view on which the 3d-plot is created
	* comp: indexes the fitness dimension for x,y and z axis in that order
	"""
	from mpl_toolkits.mplot3d import axes3d
	import matplotlib.pyplot as plt
	import numpy as np

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	fit = np.transpose([ind.cur_f for ind in pop])
	try:
		ax.plot(fit[comp[0]],fit[comp[1]],fit[comp[2]], 'ro')
	except IndexError:
		print 'Error. Please choose correct fitness dimensions for printing!'

	ax.view_init(azim=a)
	plt.show()
	return ax


def _dtlz234_plot(self, pop, a=40, comp=[0,1,2]):
	"""
	Specific plot-method for the DTLZ2, DTLZ3 and DTLZ4 - plotting also the optimal pareto-front

	USAGE: prob.plot(pop, comp[0,2,3])
	* pop: population of solutions to the problem
	* a: angle of view on which the 3d-plot is created
	* comp: indexes the fitness dimension for x,y and z axis in that order
	"""

	from mpl_toolkits.mplot3d import axes3d
	import matplotlib.pyplot as plt
	import numpy as np

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	# plot the wireframe of the known optimal pareto front
	thetas = np.linspace(0, (np.pi / 2.0), 30)
	#gammas = np.linspace(-np.pi / 4, np.pi / 4, 30)
	gammas = np.linspace(0, (np.pi / 2.0), 30)

	x_frame = np.outer(np.cos(thetas), np.cos(gammas))
	y_frame = np.outer(np.cos(thetas), np.sin(gammas))
	z_frame = np.outer(np.sin(thetas), np.ones(np.size(gammas)))

	ax.view_init(azim=a)

	ax.set_autoscalex_on(False)
	ax.set_autoscaley_on(False)
	ax.set_autoscalez_on(False)

	ax.set_xlim(0, 1.8)
	ax.set_ylim(0, 1.8)
	ax.set_zlim(0, 1.8)

	ax.plot_wireframe(x_frame,y_frame,z_frame)

	# plot the individuals of the population
	fit = np.transpose([ind.cur_f for ind in pop])
	try:
		ax.plot(fit[comp[0]],fit[comp[1]],fit[comp[2]], 'ro')
	except IndexError:
		print 'Error. Please choose correct fitness dimensions for printing!'
	plt.show()
	return ax


def _dtlz1_ctor(self, k = 5, fdim = 3):
	"""
	Constructs a DTLZ1 problem (Box-Constrained, continuous, multimodal, scalable multi-objective)

	NOTE: K Deb, L Thiele, M Laumanns, E Zitzler, Scalable test problems for evolutionary multiobjective optimization

	The optimal pareto front lies on a linear hyperplane \sum_{m=1}^M f_m = 0.5.

	USAGE: problem.dt1z1(k = 5, fdim=3)
	* k: defines problemdimension as k + fdim - 1
	* fdim: number of objectives
	"""

	arg_list=[k, fdim]
	self._orig_init(*arg_list)


dtlz1._orig_init = dtlz1.__init__
dtlz1.__init__ = _dtlz1_ctor
dtlz1.plot = _mo3d_plot


def _dtlz1_p_distance(self, pop):
	"""
	This metric is the value of the distance function g inherent to this particular multi-objective problem	averaged over the whole population. It is 0.0 iff all individuals of the population are pareto-optimal.
		
	USAGE: x = prob.p_distance(isl.population)
	
	* pop: population to evaluate
	"""
	return self._orig_p_distance(pop)

dtlz1._orig_p_distance = dtlz1.p_distance
dtlz1.p_distance = _dtlz1_p_distance


def _dtlz2_ctor(self, k = 10, fdim = 3):
	"""
	Constructs a DTLZ2 problem (Box-Constrained, continuous, unimodal, scalable multi-objective)

	NOTE: K Deb, L Thiele, M Laumanns, E Zitzler, Scalable test problems for evolutionary multiobjective optimization

	The search space is continous, unimodal and the problem is not deceptive. A rather easy problem.

	USAGE: problem.dt1z2(k = 10, fdim=3)
	* k: defines problemdimension as k + fdim - 1
	* fdim: number of objectives
	"""

	arg_list=[k, fdim]
	self._orig_init(*arg_list)

dtlz2._orig_init = dtlz2.__init__
dtlz2.__init__ = _dtlz2_ctor
dtlz2.plot = _dtlz234_plot


def _dtlz3_ctor(self, k = 10, fdim = 3):
	"""
	Constructs a DTLZ3 problem (Box-Constrained, continuous, multimodal, scalable multi-objective)

	It is supposed to be harder to converge towards the optimal pareto front than DTLZ2

	NOTE: K Deb, L Thiele, M Laumanns, E Zitzler, Scalable test problems for evolutionary multiobjective optimization

	USAGE: problem.dt1z3(k = 10, fdim=3)
	* k: defines problemdimension as k + fdim - 1
	* fdim: number of objectives
	"""

	arg_list=[k, fdim]
	self._orig_init(*arg_list)

dtlz3._orig_init = dtlz3.__init__
dtlz3.__init__ = _dtlz3_ctor
dtlz3.plot = _dtlz234_plot


def _dtlz4_ctor(self, k = 10, fdim = 3, alpha=100):
	"""
	Constructs a DTLZ4 problem (Box-Constrained, continuous, unimodal, scalable multi-objective)

	The search space contains a dense area of solutions next to the f_M/f_1 plane.

	NOTE: K Deb, L Thiele, M Laumanns, E Zitzler, Scalable test problems for evolutionary multiobjective optimization

	USAGE: problem.dt1z4(k = 10, fdim=3, alpha=100)
	* k: defines problem dimension as k + fdim - 1
	* fdim: number of objectives
	* alpha: parameter controlling density of solutions
	"""

	arg_list=[k, fdim, alpha]
	self._orig_init(*arg_list)

dtlz4._orig_init = dtlz4.__init__
dtlz4.__init__ = _dtlz4_ctor
dtlz4.plot = _dtlz234_plot


def _dtlz5_ctor(self, k = 10, fdim = 3):
	"""
	Constructs a DTLZ5 problem (Box-Constrained, continuous, unimodal, scalable multi-objective)

	This problem will test an MOEA's ability to converge to a cruve and will also allow an easier way to visually demonstrate
	(just by plotting f_M with any other objective function) the performance of an MOEA. Since there is a natural bias for
	solutions close to this Pareto-optimal curve, this problem may be easy for an algorithmn to solve. Because of its simplicity
	its recommended to use a higher number of objectives M \in [5, 10].

	NOTE: K Deb, L Thiele, M Laumanns, E Zitzler, Scalable test problems for evolutionary multiobjective optimization

	USAGE: problem.dt1z5(k = 10, fdim = 3)
	* k: defines problem dimension as k + fdim - 1
	* fdim: number of objectives
	"""

	arg_list=[k, fdim]
	self._orig_init(*arg_list)

dtlz5._orig_init = dtlz5.__init__
dtlz5.__init__ = _dtlz5_ctor
dtlz5.plot = _mo3d_plot


def _dtlz6_ctor(self, k = 10, fdim = 3):
	"""
	Constructs a DTLZ6 problem (Box-Constrained, continuous, unimodal, scalable multi-objective)

	A more difficult version of the DTLZ6 problem: the non-linear distance function g makes it harder to convergence
	against the pareto optimal curve.

	NOTE: K Deb, L Thiele, M Laumanns, E Zitzler, Scalable test problems for evolutionary multiobjective optimization

	USAGE: problem.dt1z6(k = 10, fdim = 3)
	* k: defines problem dimension as k + fdim - 1
	* fdim: number of objectives
	"""

	arg_list=[k, fdim]
	self._orig_init(*arg_list)

dtlz6._orig_init = dtlz6.__init__
dtlz6.__init__ = _dtlz6_ctor
dtlz6.plot = _mo3d_plot


def _dtlz7_ctor(self, k = 20, fdim = 3):
	"""
	Constructs a DTLZ7 problem (Box-Constrained, continuous, scalable multi-objective)

	This problem has disconnected Pareto-optimal regions in the search space.

	NOTE: K Deb, L Thiele, M Laumanns, E Zitzler, Scalable test problems for evolutionary multiobjective optimization

	USAGE: problem.dt1z7(k = 20, fdim = 3)
	* k: defines problem dimension as k + fdim - 1
	* fdim: number of objectives
	"""

	arg_list=[k, fdim]
	self._orig_init(*arg_list)

dtlz7._orig_init = dtlz7.__init__
dtlz7.__init__ = _dtlz7_ctor
dtlz7.plot = _mo3d_plot


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

def _kur_ctor(self,dim = 10):
	"""
	Constructs a Kursawe's study problem (Box-Constrained Continuous Multi-Objective)

	NOTE: K Deb, A Pratap, S Agarwal: A fast and elitist multiobjective genetic algorithm: NSGA-II, IEEE Transactions on, 2002

	USAGE: problem.kur()

	* dim: problem dimension
	"""
	arg_list=[]
	arg_list.append(dim)
	self._orig_init(*arg_list)

kur._orig_init = kur.__init__
kur.__init__ = _kur_ctor

def _fon_ctor(self):
	"""
	Constructs a Fonseca and Fleming's study problem (Box-Constrained Continuous Multi-Objective)

	NOTE: K Deb, A Pratap, S Agarwal: A fast and elitist multiobjective genetic algorithm: NSGA-II, IEEE Transactions on, 2002

	USAGE: problem.fon()
	"""
	arg_list=[]
	self._orig_init(*arg_list)

fon._orig_init = fon.__init__
fon.__init__ = _fon_ctor

def _pol_ctor(self):
	"""
	Constructs a Poloni's study study problem (Box-Constrained Continuous Multi-Objective)

	NOTE: K Deb, A Pratap, S Agarwal: A fast and elitist multiobjective genetic algorithm: NSGA-II, IEEE Transactions on, 2002

	USAGE: problem.pol()
	"""
	arg_list=[]
	self._orig_init(*arg_list)

pol._orig_init = pol.__init__
pol.__init__ = _pol_ctor

def _sch_ctor(self):
	"""
	Constructs a Schaffer's study problem (Box-Constrained Continuous Multi-Objective)

	NOTE: K Deb, A Pratap, S Agarwal: A fast and elitist multiobjective genetic algorithm: NSGA-II, IEEE Transactions on, 2002

	USAGE: problem.sch()
	"""
	arg_list=[]
	self._orig_init(*arg_list)

sch._orig_init = sch.__init__
sch.__init__ = _sch_ctor

def _zdt1_p_distance(self, pop):
	"""
	This metric is the value of the distance function g inherent to this particular multi-objective problem	averaged over the whole population. It is 0.0 iff all individuals of the population are pareto-optimal.
		
	USAGE: x = prob.p_distance(isl.population)
	
	* pop: population to evaluate
	"""
	return self._orig_p_distance(pop)

zdt1._orig_p_distance = zdt1.p_distance
zdt1.p_distance = _zdt1_p_distance


def _zdt1_ctor(self, dim = 30):
	"""
	Constructs a ZDT1 problem (Box-Constrained Continuous Multi-Objective)

	NOTE: K Deb, A Pratap, S Agarwal: A fast and elitist multiobjective genetic algorithm: NSGA-II, IEEE Transactions on, 2002

	USAGE: problem.zdt1(dim = 30)

	* dim: problem dimension
	"""

	arg_list=[]
	arg_list.append(dim)
	self._orig_init(*arg_list)
	
zdt1._orig_init = zdt1.__init__
zdt1.__init__ = _zdt1_ctor


def _zdt2_p_distance(self, pop):
	"""
	This metric is the value of the distance function g inherent to this particular multi-objective problem	averaged over the whole population. It is 0.0 iff all individuals of the population are pareto-optimal.
		
	USAGE: x = prob.p_distance(isl.population)
	
	* pop: population to evaluate
	"""
	return self._orig_p_distance(pop)

zdt2._orig_p_distance = zdt2.p_distance
zdt2.p_distance = _zdt2_p_distance


def _zdt2_ctor(self, dim = 30):
	"""
	Constructs a ZDT2 problem (Box-Constrained Continuous Multi-Objective)

	NOTE: K Deb, A Pratap, S Agarwal: A fast and elitist multiobjective genetic algorithm: NSGA-II, IEEE Transactions on, 2002

	USAGE: problem.zdt2(dim = 30)

	* dim: problem dimension
	"""
	arg_list=[]
	arg_list.append(dim)
	self._orig_init(*arg_list)

zdt2._orig_init = zdt2.__init__
zdt2.__init__ = _zdt2_ctor


def _zdt3_p_distance(self, pop):
	"""
	This metric is the value of the distance function g inherent to this particular multi-objective problem	averaged over the whole population. It is 0.0 iff all individuals of the population are pareto-optimal.
		
	USAGE: x = prob.p_distance(isl.population)
	
	* pop: population to evaluate
	"""
	return self._orig_p_distance(pop)

zdt3._orig_p_distance = zdt3.p_distance
zdt3.p_distance = _zdt3_p_distance


def _zdt3_ctor(self, dim = 30):
	"""
	Constructs a ZDT3 problem (Box-Constrained Continuous Multi-Objective)

	NOTE: K Deb, A Pratap, S Agarwal: A fast and elitist multiobjective genetic algorithm: NSGA-II, IEEE Transactions on, 2002

	USAGE: problem.zdt3(dim = 30)

	* dim: problem dimension
	"""
	arg_list=[]
	arg_list.append(dim)
	self._orig_init(*arg_list)

zdt3._orig_init = zdt3.__init__
zdt3.__init__ = _zdt3_ctor


def _zdt4_p_distance(self, pop):
	"""
	This metric is the value of the distance function g inherent to this particular multi-objective problem	averaged over the whole population. It is 0.0 iff all individuals of the population are pareto-optimal.
		
	USAGE: x = prob.p_distance(isl.population)
	
	* pop: population to evaluate
	"""
	return self._orig_p_distance(pop)

zdt4._orig_p_distance = zdt4.p_distance
zdt4.p_distance = _zdt4_p_distance


def _zdt4_ctor(self, dim = 10):
	"""
	Constructs a ZDT4 problem (Box-Constrained Continuous Multi-Objective)

	NOTE: K Deb, A Pratap, S Agarwal: A fast and elitist multiobjective genetic algorithm: NSGA-II, IEEE Transactions on, 2002

	USAGE: problem.zdt4(dim = 10)

	* dim: problem dimension
	"""
	arg_list=[]
	arg_list.append(dim)
	self._orig_init(*arg_list)

zdt4._orig_init = zdt4.__init__
zdt4.__init__ = _zdt4_ctor


def _zdt6_p_distance(self, pop):
	"""
	This metric is the value of the distance function g inherent to this particular multi-objective problem	averaged over the whole population. It is 0.0 iff all individuals of the population are pareto-optimal.
		
	USAGE: x = prob.p_distance(isl.population)
	
	* pop: population to evaluate
	"""
	return self._orig_p_distance(pop)

zdt6._orig_p_distance = zdt6.p_distance
zdt6.p_distance = _zdt6_p_distance


def _zdt6_ctor(self, dim = 10):
	"""
	Constructs a ZDT6 problem (Box-Constrained Continuous Multi-Objective)

	NOTE: K Deb, A Pratap, S Agarwal: A fast and elitist multiobjective genetic algorithm: NSGA-II, IEEE Transactions on, 2002

	USAGE: problem.zdt6(dim = 10)

	* dim: problem dimension
	"""
	arg_list=[]
	arg_list.append(dim)
	self._orig_init(*arg_list)

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
	self._orig_init(*arg_list)

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


