from _astro_toolbox import *

def lambertI(r0,r1,T,mu,way):
	from PyGMO import vector
	from _astro_toolbox import __lambertI, __lambert_result
	lr = __lambert_result()
	__lambertI(r0,r1,T,mu,way,lr)
	return lr.v0, lr.v1, lr.a, lr.p, lr.theta, lr.it
