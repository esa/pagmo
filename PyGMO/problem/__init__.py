# -*- coding: iso-8859-1 -*-
from _problem import *

# Raw C++ base class.
_base = _problem._base

class base(_problem._base):
	def __init__(self,*args):
		if len(args) == 0:
			raise ValueError("Cannot initialise base problem without parameters for the constructor.")
		_problem._base.__init__(self,*args)
	def _get_typename(self):
		return str(type(self))
	def __get_deepcopy__(self):
		from copy import deepcopy
		return deepcopy(self)
	def get_name(self):
		return self._get_typename()

class py_test(base):
	"""
	Minimal test problem implemented purely in Python. Objective function
	is De Jong's unidimensional sphere function.
	"""
	def __init__(self):
		super(py_test,self).__init__(1)
	def _objfun_impl(self,x):
		return (x[0] * x[0],)
		


class py_earth_mars(base):
	"""
	This is a PaGMO.problem object that represents a low-thrust transfer between Earth and Mars.
	"""
	def __init__(self,mass=1000,Tmax=0.05,Isp=2500,Vinf=3,nseg=10):
		"""The constructor first calls the base class constructor and assignes values to data members"""
		super(py_earth_mars,self).__init__(6 + nseg*3,0,1,8 + nseg,nseg+1,1e-7)
		try:
			import PyKEP
			self.__earth = PyKEP.planet_ss('earth')
			self.__mars = PyKEP.planet_ss('mars')
			self.__sc = PyKEP.sims_flanagan.spacecraft(mass,Tmax,Isp)
			self.__Vinf = Vinf*1000
			self.__leg = PyKEP.sims_flanagan.leg()
			self.__leg.set_mu(PyKEP.MU_SUN)
			self.__leg.set_sc(self.__sc)
			self.__nseg = nseg
			self.set_bounds([0,60,self.__sc.mass/10,-self.__Vinf,-self.__Vinf,-self.__Vinf] + [-1] * 3 *nseg,[3000,1500,self.__sc.mass,self.__Vinf,self.__Vinf,self.__Vinf] + [1] * 3 * nseg)
		except ImportError:
			raise ImportError('PyKEP is needed to instantiate an earth_mass')
	def _objfun_impl(self,x):
		return (-x[2],)
	def _compute_constraints_impl(self,x):
		import PyKEP
		start = PyKEP.epoch(x[0])
		end = PyKEP.epoch(x[0] + x[1])
		r,v = self.__earth.eph(start)
		v_list = list(v)
		v_list[0] += x[3]
		v_list[1] += x[4]
		v_list[2] += x[5]
		x0 = PyKEP.sims_flanagan.sc_state(r,v_list,self.__sc.mass)
		r,v = self.__mars.eph(end)
		xe = PyKEP.sims_flanagan.sc_state(r, v ,x[2])
		self.__leg.set(start,x0,x[-3 * self.__nseg:],end,xe)
		v_inf_con = (x[3] * x[3] + x[4] * x[4] + x[5] * x[5] - self.__Vinf * self.__Vinf) / (PyKEP.EARTH_VELOCITY * PyKEP.EARTH_VELOCITY)
		retval = list(self.__leg.mismatch_constraints() + self.__leg.throttles_constraints()) + [v_inf_con]
		retval[0] /= PyKEP.AU
		retval[1] /= PyKEP.AU
		retval[2] /= PyKEP.AU
		retval[3] /= PyKEP.EARTH_VELOCITY
		retval[4] /= PyKEP.EARTH_VELOCITY
		retval[5] /= PyKEP.EARTH_VELOCITY
		retval[6] /= self.__sc.mass
		return retval
	def high_fidelity(self,state):
		self.__leg.high_fidelity(state)



#class apophis_impact(base):
	#def __init__(self):
		#super(apophis_impact,self).__init__(2)
		#import PyKEP
		#import numpy
		#self.set_bounds([PyKEP.epoch_from_string('2020-06-01 00:00:00').mjd2000(),30],
			#[PyKEP.epoch_from_string('2021-06-01 00:00:00').mjd2000(),1000])
	#def _objfun_impl(self,x):
		#import PyKEP
		#from numpy import array, dot
		#from numpy.linalg import norm
		#from math import exp
		#earth = PyKEP.planet_ss('earth')
		#apophis = PyKEP.planet_mpcorb('99942   19.2   0.15 K107N 202.49545  126.41859  204.43202    3.33173  0.1911104  1.11267324   0.9223398  1 MPO164109  1397   2 2004-2008 0.40 M-v 3Eh MPCAPO     C802  (99942) Apophis            20080109')
		#r0, v0 = earth.eph(PyKEP.epoch(x[0]))
		#rf, vf = apophis.eph(PyKEP.epoch(x[0] + x[1]))
		#leg = PyKEP.lambert_problem(r0,rf,x[1] * PyKEP.DAY2SEC,PyKEP.MU_SUN)
		#vsc_1 = array(leg.get_v1()[0])
		#vsc_2 = array(leg.get_v2()[0])
		#U = array(vf) - vsc_2
		#vU = dot(array(vf),U)
		#DV0 = max(0.,norm(vsc_1 - array(v0)) - 4472.1359549995796)
		##DV0 = max(0.,norm(vsc_1 - array(v0)) - 3500.)
		#m = 500. * exp(-DV0 / 233. / PyKEP.G0)
		#return (-abs(m * vU),)
	#def pretty(self,x):
		#import PyKEP
		#from numpy import array, dot
		#from numpy.linalg import norm
		#from math import exp
		#earth = PyKEP.planet_ss('earth')
		#apophis = PyKEP.planet_mpcorb('99942   19.2   0.15 K107N 202.49545  126.41859  204.43202    3.33173  0.1911104  1.11267324   0.9223398  1 MPO164109  1397   2 2004-2008 0.40 M-v 3Eh MPCAPO     C802  (99942) Apophis            20080109')
		#r0, v0 = earth.eph(PyKEP.epoch(x[0]))
		#rf, vf = apophis.eph(PyKEP.epoch(x[0] + x[1]))
		#leg = PyKEP.lambert_problem(r0,rf,x[1] * PyKEP.DAY2SEC,PyKEP.MU_SUN)
		#vsc_1 = array(leg.get_v1()[0])
		#vsc_2 = array(leg.get_v2()[0])
		#U = array(vf) - vsc_2
		#vU = dot(array(vf),U)
		#DV0 = max(0.,norm(vsc_1 - array(v0)) - 4472.1359549995796)
		##DV0 = max(0.,norm(vsc_1 - array(v0)) - 3500.)
		#m = 500. * exp(-DV0 / 233. / PyKEP.G0)
		#print('Earth eph: ' + str(r0) + ',' + str(v0))
		#print('Apophis eph: ' + str(rf) + ',' + str(vf))
		#print('SC velocities: ' + str(vsc_1) + ',' + str(vsc_2))
		#print('Rel. velocity at impact: ' + str(U))
		#print('DV: ' + str(DV0))
		#print('Final mass: ' + str(m))
		#print('Starting date: ' + str(PyKEP.epoch(x[0])))
		#print('Arrival date: ' + str(PyKEP.epoch(x[0] + x[1])))

def _get_problem_list():
	from PyGMO import problem
	return [problem.__dict__[n] for n in filter(lambda n: not n.startswith('_') and not n == 'base' and issubclass(problem.__dict__[n],problem._base),dir(problem))]
