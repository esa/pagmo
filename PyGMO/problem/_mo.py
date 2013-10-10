from _problem import zdt, zdt1, zdt2, zdt3, zdt4, zdt5, zdt6, dtlz, dtlz1, dtlz2, dtlz3, dtlz4, dtlz5, dtlz6, dtlz7

def _mo3d_plot(pop, a=40, comp=[0,1,2]):
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




def _dtlz234_plot(pop, a=40, comp=[0,1,2]):
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

def _zdt_ctor(self, prob_id = 1, param_1 = None):
        """
        Constructs a multi-objective box-constrained problem from the ZDT testsuite
        
        NOTE: K Deb, A Pratap, S Agarwal: A fast and elitist multiobjective genetic algorithm: NSGA-II, IEEE Transactions on, 2002
        
        USAGE: problem.zdt(prob_id = 1, param_1 = 30)
        
        * prob_id: Problem number, one of [1,2,...24]
        * param_1: problem dimension for all ZDT problems except ZDT5 (here it is the number of binary strings used)
        """

	arg_list=[]
	arg_list.append(prob_id)
	if param_1 == None:
		if prob_id in [1,2,3,4]:
			arg_list.append(30)
		elif prob_id == 5:
			arg_list.append(11)
		else:
			arg_list.append(10)
	else:
		arg_list.append(param_1)
        self._orig_init(*arg_list)
        
zdt._orig_init = zdt.__init__
zdt.__init__ = _zdt_ctor



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

def _zdt5_ctor(self, bnr = 11):
	"""
	Constructs a ZDT5 problem (Box-Constrained Continuous Multi-Objective)

	NOTE: E. Zitzler, K. Deb, L. Thiele: Comparison of Multiobjective Evolutionary Algorithms: Empirical Results, Evolutionary Computation, 2000

	USAGE: problem.zdt5(bnr = 11)

	* bnr: number of binary strings used for the problem (problem dimension is 30 + 5 * (bnr - 1))
	"""
	arg_list=[]
	arg_list.append(bnr)
	self._orig_init(*arg_list)

zdt5._orig_init = zdt5.__init__
zdt5.__init__ = _zdt5_ctor

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

def _dtlz_ctor(self, prob_id = 1, k = None, fdim = 3, alpha = 100):
        """
        Constructs a multi-objective box-constrained problem from the DTLZ testsuite
        
	NOTE: K Deb, L Thiele, M Laumanns, E Zitzler, Scalable test problems for evolutionary multiobjective optimization
        
        USAGE: problem.dtlz(prob_id = 1, k = 20, fdim = 4)
        
        * prob_id: Problem number, one of [1,2,...7]
        * k: paramter defining integer dimension of the problem: k + fdim - 1
        * fdim: number of objectives
        * alpha: controls density of solutions (just used for prob_id = 4)
        """

	arg_list=[]
	arg_list.append(prob_id)
	if k == None:
		if prob_id == 1:
			arg_list.append(5)
		elif prob_id in [2,3,4,5,6]:
			arg_list.append(10)
		else:
			arg_list.append(20)
	else:
		arg_list.append(k)
	arg_list.append(fdim)
	arg_list.append(alpha)
        self._orig_init(*arg_list)
        if prob_id in [2,3,4]:
	    self.plot = _dtlz234_plot
	else:
	    self.plot = _mo3d_plot
        
dtlz._orig_init = dtlz.__init__
dtlz.__init__ = _dtlz_ctor


def _dtlz1_ctor(self, k = 5, fdim = 3):
	"""
	Constructs a DTLZ1 problem (Box-Constrained, continuous, multimodal, scalable multi-objective)

	NOTE: K Deb, L Thiele, M Laumanns, E Zitzler, Scalable test problems for evolutionary multiobjective optimization

	The optimal pareto front lies on a linear hyperplane with f_1 + f_2 + ... + f_m = 0.5.

	USAGE: problem.dt1z1(k = 5, fdim=3)
	
	* k: defines problem dimension as k + fdim - 1
	
	* fdim: number of objectives
	"""

	arg_list=[k, fdim]
	self._orig_init(*arg_list)
	self.plot = _mo3d_plot

dtlz1._orig_init = dtlz1.__init__
dtlz1.__init__ = _dtlz1_ctor


def _dtlz2_ctor(self, k = 10, fdim = 3):
	"""
	Constructs a DTLZ2 problem (Box-Constrained, continuous, unimodal, scalable multi-objective)

	NOTE: K Deb, L Thiele, M Laumanns, E Zitzler, Scalable test problems for evolutionary multiobjective optimization

	The search space is continous, unimodal and the problem is not deceptive. The Pareto-front is a quarter of a sphere

	USAGE: problem.dt1z2(k = 10, fdim=3)
	
	* k: defines problemdimension as k + fdim - 1
	
	* fdim: number of objectives
	"""

	arg_list=[k, fdim]
	self._orig_init(*arg_list)
	self.plot = _dtlz234_plot

dtlz2._orig_init = dtlz2.__init__
dtlz2.__init__ = _dtlz2_ctor


def _dtlz3_ctor(self, k = 10, fdim = 3):
	"""
	Constructs a DTLZ3 problem (Box-Constrained, continuous, multimodal, scalable multi-objective)

	Same Pareto Front as DTLZ2, but harder to converge towards it

	NOTE: K Deb, L Thiele, M Laumanns, E Zitzler, Scalable test problems for evolutionary multiobjective optimization

	USAGE: problem.dt1z3(k = 10, fdim=3)
	
	* k: defines problemdimension as k + fdim - 1
	
	* fdim: number of objectives
	"""

	arg_list=[k, fdim]
	self._orig_init(*arg_list)
	self.plot = _dtlz234_plot

dtlz3._orig_init = dtlz3.__init__
dtlz3.__init__ = _dtlz3_ctor


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
	self.plot = _dtlz234_plot

dtlz4._orig_init = dtlz4.__init__
dtlz4.__init__ = _dtlz4_ctor


def _dtlz5_ctor(self, k = 10, fdim = 3):
	"""
	Constructs a DTLZ5 problem (Box-Constrained, continuous, unimodal, scalable multi-objective)

	This problem will test an MOEA's ability to converge to a curve and will also allow an easier way to visually demonstrate
	(just by plotting f_M with any other objective function) the performance of an MOEA. Since there is a natural bias for
	solutions close to this Pareto-optimal curve, this problem may be easy for an algorithmn to solve. Because of its simplicity
	its recommended to use a higher number of objectives.

	NOTE: K Deb, L Thiele, M Laumanns, E Zitzler, Scalable test problems for evolutionary multiobjective optimization

	USAGE: problem.dt1z5(k = 10, fdim = 3)
	
	* k: defines problem dimension as k + fdim - 1
	
	* fdim: number of objectives
	"""

	arg_list=[k, fdim]
	self._orig_init(*arg_list)
	self.plot = _mo3d_plot

dtlz5._orig_init = dtlz5.__init__
dtlz5.__init__ = _dtlz5_ctor

def _dtlz6_ctor(self, k = 10, fdim = 3):
	"""
	Constructs a DTLZ6 problem (Box-Constrained, continuous, unimodal, scalable multi-objective)

	A more difficult version of the DTLZ5 problem: the non-linear distance function g makes it harder to convergence
	against the pareto optimal curve.

	NOTE: K Deb, L Thiele, M Laumanns, E Zitzler, Scalable test problems for evolutionary multiobjective optimization

	USAGE: problem.dt1z6(k = 10, fdim = 3)
	
	* k: defines problem dimension as k + fdim - 1
	
	* fdim: number of objectives
	"""

	arg_list=[k, fdim]
	self._orig_init(*arg_list)
	self.plot = _mo3d_plot

dtlz6._orig_init = dtlz6.__init__
dtlz6.__init__ = _dtlz6_ctor

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
	self.plot = _mo3d_plot

dtlz7._orig_init = dtlz7.__init__
dtlz7.__init__ = _dtlz7_ctor