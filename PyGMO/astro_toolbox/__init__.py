from _astro_toolbox import *

def lambertI(r0,r1,T,mu,way):
	"""
	Lambert.
	"""
	from _astro_toolbox import __lambertI, __lambert_result
	lr = __lambert_result()
	__lambertI(r0,r1,T,mu,way,lr)
	return lr.v0, lr.v1, lr.a, lr.p, lr.theta, lr.it

def propagate_kep(r0,v0,t,mu):
	"""
	Keplerian propagation.
	"""
	from PyGMO import vector
	from _astro_toolbox import __propagate_kep
	r1 = vector([0.,0.,0.])
	v1 = vector([0.,0.,0.])
	__propagate_kep(r0,v0,t,mu,r1,v1)
	return r1, v1
