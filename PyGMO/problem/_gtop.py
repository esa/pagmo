from _problem import *
from _problem import _gtoc_2_objective 
		
# Redefining the constructors of all problems to obtain good documentation and allowing kwargs
def _cassini_1_ctor(self):
	"""
	Constructs a Cassini 1 Problem (Box-Constrained Continuous Single-Objective)
	
	NOTE: This problem (MGA) belongs to the GTOP database [http://www.esa.int/gsp/ACT/inf/op/globopt.htm]
	      It has a global minimum at 4.9307 [km/s], and it is a deceptive problem
	      with a larger minimum at 5.303 [km/s]
	      
	USAGE: problem.cassini_1()
	
	"""
	
	# We construct the arg list for the original constructor exposed by boost_python
	arg_list=[]
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
