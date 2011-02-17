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


class py_pl2pl(base):
	"""
	This problem represents a low-thrust transfer between a departure planet (default is the earth)
	and a target planet (default is mars). The spacecraft is described 
	by its starting mass (mass) its engine specific impulse (Isp) and its engine maximum thrust (Tmax). The
	Sims-Flanagan model is used to describe a trajectory. A variable number of segments (nseg) can be used
	An initial velocity with respect to the Earth is allowed (Vinf) assumed to be given by the launcher
	The method high_fidelity allows to use a continuous thrust model rather than impulses
	"""
	try:
		import PyKEP
		earth = PyKEP.planet_ss('earth')
		mars = PyKEP.planet_ss('mars')
		def __init__(self,mass=1000,Tmax=0.05,Isp=2500,Vinf_0=3,Vinf_f=0,nseg=10,departure = earth, target = mars):
			"""__init__(self,mass=1000,Tmax=0.05,Isp=2500,Vinf=3,nseg=10,departure = earth, target = mars)"""
			import PyKEP
			super(py_pl2pl,self).__init__(9 + nseg*3,0,1,9 + nseg,nseg+2,1e-5)
			self.__departure = departure
			self.__target = target
			self.__sc = PyKEP.sims_flanagan.spacecraft(mass,Tmax,Isp)
			self.__Vinf_0 = Vinf_0*1000
			self.__Vinf_f = Vinf_f*1000
			self.__leg = PyKEP.sims_flanagan.leg()
			self.__leg.set_mu(PyKEP.MU_SUN)
			self.__leg.set_spacecraft(self.__sc)
			self.__nseg = nseg
			self.set_bounds([0,60,self.__sc.mass/10,-self.__Vinf_0,-self.__Vinf_0,-self.__Vinf_0,-self.__Vinf_f,-self.__Vinf_f,-self.__Vinf_f] + [-1] * 3 *nseg,[3000,1500,self.__sc.mass,self.__Vinf_0,self.__Vinf_0,self.__Vinf_0,self.__Vinf_f,self.__Vinf_f,self.__Vinf_f] + [1] * 3 * nseg)
		def _objfun_impl(self,x):
			return (-x[2],)
		def _compute_constraints_impl(self,x):
			import PyKEP
			start = PyKEP.epoch(x[0])
			end = PyKEP.epoch(x[0] + x[1])
			r,v = self.__departure.eph(start)
			v_list = list(v)
			v_list[0] += x[3]
			v_list[1] += x[4]
			v_list[2] += x[5]
			x0 = PyKEP.sims_flanagan.sc_state(r,v_list,self.__sc.mass)
			r,v = self.__target.eph(end)
			xe = PyKEP.sims_flanagan.sc_state(r, v ,x[2])
			self.__leg.set(start,x0,x[-3 * self.__nseg:],end,xe)
			v_inf_con_0 = (x[3] * x[3] + x[4] * x[4] + x[5] * x[5] - self.__Vinf_0 * self.__Vinf_0) / (PyKEP.EARTH_VELOCITY * PyKEP.EARTH_VELOCITY)
			v_inf_con_f = (x[6] * x[6] + x[7] * x[7] + x[8] * x[8] - self.__Vinf_f * self.__Vinf_f) / (PyKEP.EARTH_VELOCITY * PyKEP.EARTH_VELOCITY)
			retval = list(self.__leg.mismatch_constraints() + self.__leg.throttles_constraints()) + [v_inf_con_0] + [v_inf_con_f] 
			retval[0] /= PyKEP.AU
			retval[1] /= PyKEP.AU
			retval[2] /= PyKEP.AU
			retval[3] /= PyKEP.EARTH_VELOCITY
			retval[4] /= PyKEP.EARTH_VELOCITY
			retval[5] /= PyKEP.EARTH_VELOCITY
			retval[6] /= self.__sc.mass
			return retval
		def pretty(self,x):
			"""Decodes the decision vector x"""
			import PyKEP
			start = PyKEP.epoch(x[0])
			end = PyKEP.epoch(x[0] + x[1])
			r,v = self.__departure.eph(start)
			v_list = list(v)
			v_list[0] += x[3]
			v_list[1] += x[4]
			v_list[2] += x[5]
			x0 = PyKEP.sims_flanagan.sc_state(r,v_list,self.__sc.mass)
			r,v = self.__target.eph(end)
			xe = PyKEP.sims_flanagan.sc_state(r, v ,x[2])
			self.__leg.set(start,x0,x[-3 * self.__nseg:],end,xe)
			print("A direct interplantary transfer\n")
			print("FROM:")
			print(self.__departure)
			print("TO:")
			print(self.__target)
			print(self.__leg)
		def get_hf(self):
			return self.__leg.high_fidelity
		def set_hf(self,state):
			self.__leg.high_fidelity = state
		high_fidelity = property(get_hf,set_hf)
	except:
		print('warning: PyKEP is not detected, skipping all PyKEP related imports')


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
