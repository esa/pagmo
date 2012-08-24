from _problem import cassini_1, gtoc_1,cassini_2, rosetta, messenger_full, tandem, laplace, sagas, gtoc_2, mga_1dsm
from _problem import _gtoc_2_objective 
		
# Redefining the constructors of all problems to obtain good documentation and allowing kwargs
def _cassini_1_ctor(self, objectives = 1):
	"""
	Constructs a Cassini 1 Problem (Box-Constrained Continuous Single-Objective)
	
	NOTE: This problem (MGA) belongs to the GTOP database [http://www.esa.int/gsp/ACT/inf/op/globopt.htm]
	      Its single objective version has a global minimum at 4.9307 [km/s],
	      and it is a deceptive problem with a larger minimum at 5.303 [km/s]
	      
	USAGE: problem.cassini_1(objectives = 1)
	
	* objectives: number of objectives. 1=DV, 2=DV,DT
	
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(objectives)
	self._orig_init(*arg_list)
cassini_1._orig_init = cassini_1.__init__
cassini_1.__init__ = _cassini_1_ctor

def _gtoc_1_ctor(self):
	"""
	Constructs a GTOC 1 Problem (Box-Constrained Continuous Single-Objective)
	
	NOTE: This problem (MGA) belongs to the GTOP database [http://www.esa.int/gsp/ACT/inf/op/globopt.htm]
	      
	      Best known global minima is at -1,581,950
	      
	USAGE: problem.gtoc_1()
	
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	self._orig_init(*arg_list)
gtoc_1._orig_init = gtoc_1.__init__
gtoc_1.__init__ = _gtoc_1_ctor

def _cassini_2_ctor(self):
	"""
	Constructs a Cassini 2 Problem (Box-Constrained Continuous Single-Objective)
	
	NOTE: This problem (MGA-1DSM) belongs to the GTOP database [http://www.esa.int/gsp/ACT/inf/op/globopt.htm]
	      It models the same interplanetary trajectory as the cassini_1 problem, but
	      in a more accurate fashion, allowing deep space manouvres
	      
	      Best known global minimum is at 8.383 [km/s]
	      
	USAGE: problem.cassini_2()
	
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	self._orig_init(*arg_list)
cassini_2._orig_init = cassini_2.__init__
cassini_2.__init__ = _cassini_2_ctor

def _rosetta_ctor(self):
	"""
	Constructs a Rosetta Problem (Box-Constrained Continuous Single-Objective)
	
	NOTE: This problem (MGA-1DSM) belongs to the GTOP database [http://www.esa.int/gsp/ACT/inf/op/globopt.htm]
	      
	      Best known global minimum is at 1.343 [km/s]
	      
	USAGE: problem.rosetta()
	
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	self._orig_init(*arg_list)
rosetta._orig_init = rosetta.__init__
rosetta.__init__ = _rosetta_ctor

def _messenger_full_ctor(self):
	"""
	Constructs a Mesenger Full Problem (Box-Constrained Continuous Single-Objective)
	
	NOTE: This problem (MGA-1DSM) belongs to the GTOP database [http://www.esa.int/gsp/ACT/inf/op/globopt.htm]
	      
	      Best known global minimum is at 2.970 [km/s], but physics indicate a minimum exist at 2.3 ....
	      
	USAGE: problem.messenger_full()
	
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	self._orig_init(*arg_list)
messenger_full._orig_init = messenger_full.__init__
messenger_full.__init__ = _messenger_full_ctor

def _tandem_ctor(self, prob_id = 7, max_tof = -1):
	"""
	Constructs a TandEM Problem (Box-Constrained Continuous Single-Objective)
	
	NOTE: This problem (MGA-1DSM) belongs to the GTOP database [http://www.esa.int/gsp/ACT/inf/op/globopt.htm]. The objective function is -log(m_final).
	      
	USAGE: problem.tandem(prob_id = 7, max_tof = -1)
	
	* prob_id: Selects the problem variant (one of 1..25). All problems differ from the fly-by sequence
	* max_tof = Activates a constriants on the maximum time of flight allowed (in years)
	
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(prob_id)
	arg_list.append(max_tof)
	self._orig_init(*arg_list)
tandem._orig_init = tandem.__init__
tandem.__init__ = _tandem_ctor

def _laplace_ctor(self, seq = [3,2,3,3,5]):
	"""
	Constructs a EJSM-Laplace Problem (Box-Constrained Continuous Single-Objective)
	
	NOTE: This problem (MGA-1DSM) is similar to TandEM, but targets Jupiter and the user
	      can specify explicitly the planetary fly-by sequence
	      
	USAGE: problem.laplace(seq = [3,2,3,3,5])
	
	* seq: The planetary sequence. This is a list of ints that represent the planets to visit
	       1 - Mercury, 2 - Venus, 3 - Earth, 4 - Mars, 5 - Jupiter, 6 - Saturn. It must start from 3 (Earth)
	       and end with 5 (Jupiter)
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(seq)
	self._orig_init(*arg_list)
laplace._orig_init = laplace.__init__
laplace.__init__ = _laplace_ctor

def _sagas_ctor(self):
	"""
	Constructs a SAGAS Problem (Box-Constrained Continuous Single-Objective)
	
	NOTE: This problem (MGA-1DSM) belongs to the GTOP database [http://www.esa.int/gsp/ACT/inf/op/globopt.htm]
	      
	USAGE: problem.sagas()
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(seq)
	self._orig_init(*arg_list)
sagas._orig_init = sagas.__init__
sagas.__init__ = _sagas_ctor

gtoc_2.obj = _gtoc_2_objective
def _gtoc_2_ctor(self, ast1 = 815, ast2 = 300, ast3 = 110, ast4 = 47, n_seg = 10,  objective = gtoc_2.obj.MASS_TIME):
	"""
	Constructs a GTOC 2 Problem (Constrained Continuous Single-Objective)

	NOTE: This problem is a quite faitful transcription of the problem used during the GTOC2 competition
	      It Transcribe the whole OCP resulting from the low-thrust dynamics into an NLP. As such it is very
	      difficult to find feasible solutions. Note that by default the asteroid sequence isa the winning one
	      from Turin University.

	USAGE: problem.gtoc_2(ast1 = 815, ast2 = 300, ast3 = 110, ast4 = 47, n_seg = 10, objective = gtoc_2.obj.MASS_TIME)

	* ast1 id of the first asteroid to visit (Group 1:   0 - 95)
	* ast2 id of the second asteroid to visit (Group 2:  96 - 271)
	* ast3 id of the third asteroid to visit (Group 3: 272 - 571)
	* ast4 id of the fourth asteroid to visit (Group 4: 572 - 909)
	* n_seg number of segments to be used per leg
	* obj objective function in the enum {MASS,TIME,MASS_TIME}
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
	arg_list.append(ast1)
	arg_list.append(ast2)
	arg_list.append(ast3)
	arg_list.append(ast4)
	arg_list.append(n_seg)
	arg_list.append(objective)
	self._orig_init(*arg_list)
gtoc_2._orig_init = gtoc_2.__init__
gtoc_2.__init__ = _gtoc_2_ctor

try:
	from PyKEP import __version__ as d
	del d;
	#Plot of the trajectory for an mga_1dsm problem
	def _mga_1dsm_plot(self,x):
		"""
		Plots the trajectory represented by the decision vector x
		"""
		import matplotlib as mpl
		from mpl_toolkits.mplot3d import Axes3D
		import matplotlib.pyplot as plt
		from PyKEP.orbit_plots import plot_planet, plot_lambert, plot_kepler
		from PyKEP import epoch, propagate_lagrangian, lambert_problem,fb_prop, AU, MU_SUN, DAY2SEC
		from math import pi, acos, cos, sin
		from scipy.linalg import norm

		mpl.rcParams['legend.fontsize'] = 10
		fig = plt.figure()
		ax = fig.gca(projection='3d')
		ax.scatter(0,0,0, color='y')
		
		seq = self.get_sequence()
		
		n = len(x)/4
		#1 -  we 'decode' the chromosome recording the various times of flight (days) in the list T
		T = list([0]*(n))
		#a[-i] = x[-1-(i-1)*4]
		for i in xrange(n-1):	
			j = i+1;
			T[-j] = (x[5] - sum(T[-(j-1):])) * x[-1-(j-1)*4]
		T[0] = x[5] - sum(T)
		
		#2 - We compute the epochs and ephemerides of the planetary encounters
		t_P = list([None] * (n+1))
		r_P = list([None] * (n+1))
		v_P = list([None] * (n+1))
		DV = list([None] * (n+1))
		
		for i,planet in enumerate(seq):
			t_P[i] = epoch(x[0] + sum(T[0:i]))
			r_P[i],v_P[i] = planet.eph(t_P[i])
			plot_planet(ax, planet, t0=t_P[i], color=(0.8,0.6,0.8), legend=True, units = AU)

		#3 - We start with the first leg
		theta = 2*pi*x[1]
		phi = acos(2*x[2]-1)-pi/2

		Vinfx = x[3]*cos(phi)*cos(theta)
		Vinfy =	x[3]*cos(phi)*sin(theta)
		Vinfz = x[3]*sin(phi)

		v0 = [a+b for a,b in zip(v_P[0],[Vinfx,Vinfy,Vinfz])]
		r,v = propagate_lagrangian(r_P[0],v0,x[4]*T[0]*DAY2SEC,MU_SUN)
		plot_kepler(ax,r_P[0],v0,x[4]*T[0]*DAY2SEC,MU_SUN,N = 100, color='b', legend=False, units = AU)

		#Lambert arc to reach seq[1]
		dt = (1-x[4])*T[0]*DAY2SEC
		l = lambert_problem(r,r_P[1],dt,MU_SUN)
		plot_lambert(ax,l, sol = 0, color='r', legend=False, units = AU)
		v_end_l = l.get_v2()[0]
		v_beg_l = l.get_v1()[0]

		#First DSM occuring at time nu1*T1
		DV[0] = norm([a-b for a,b in zip(v_beg_l,v)])

		#4 - And we proceed with each successive leg
		for i in range(1,n):
			#Fly-by 
			v_out = fb_prop(v_end_l,v_P[i],x[7+(i-1)*4]*seq[i].radius,x[6+(i-1)*4],seq[i].mu_self)
			#s/c propagation before the DSM
			r,v = propagate_lagrangian(r_P[i],v_out,x[8+(i-1)*4]*T[i]*DAY2SEC,MU_SUN)
			plot_kepler(ax,r_P[i],v_out,x[8+(i-1)*4]*T[i]*DAY2SEC,MU_SUN,N = 100, color='b', legend=False, units = AU)
			#Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
			dt = (1-x[8+(i-1)*4])*T[i]*DAY2SEC
			l = lambert_problem(r,r_P[i+1],dt,MU_SUN)
			plot_lambert(ax,l, sol = 0, color='r', legend=False, units = AU)
			v_end_l = l.get_v2()[0]
			v_beg_l = l.get_v1()[0]
			#DSM occuring at time nu2*T2
			DV[i] = norm([a-b for a,b in zip(v_beg_l,v)])

		plt.show()
	mga_1dsm.plot = _mga_1dsm_plot
	
except:
	pass
