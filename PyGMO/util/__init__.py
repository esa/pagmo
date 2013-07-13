# -*- coding: utf-8 -*-
from _util import hypervolume, hv_algorithm
from _util.hv_algorithm import lebmeasure, optimal2d, optimal3d, wfg

__all__ = ['hypervolume', 'hv_algorithm']

def _hypervolume_ctor(self, *args, **kwargs):
	"""
	Constructs a hypervolume object used for the computation of hypervolue and exclusive hypervolume.
	Constructors use unnamed arguments.

	Object can be constructed from the population object, or from a fixed list of points
	Points within a fixed list must all be of equal size, larger than 1.

	USAGE:
		from PyGMO import *
		hv = hypervolume(pop)
		hv = hypervolume([[1,1,2],[2,1,2],[2,2,3]])

	"""
	return self._original_init(*args, **kwargs)

hypervolume._original_init = hypervolume.__init__
hypervolume.__init__ = _hypervolume_ctor

def _default_best_algorithm(f_dim):
	"""
	Choose the best suited method for given dimension.
	"""
	if f_dim == 2:
		return optimal2d()
	elif f_dim == 3:
		return optimal3d()
	else:
		return wfg()

def _hypervolume_compute(self, r, algorithm = None):
	"""
	Compute the hypervolume indicator for a given reference point, using the provided hypervolume algorithm.

	USAGE:
		hv.compute(r=[5.0]*2)
		hv.compute(r=[5.0]*2, algorithm = hv_algorithm.optimal2d())
		* r - reference point
		* algorithm - optional argument: hypervolume algorithm used for the computation.
	"""
	if not isinstance(r, list):
		raise TypeError("Reference point must be a list of lists of real numbers, e.g.: r = [1.0, 1.0, 1.0]")
	try:
		r = [float(ri) for ri in r]
	except ValueError:
		raise TypeError("Every item in reference point (r), must be castable to float, e.g.: r = [1, '2.5', 10e-4]")

	if (algorithm == None):
		algorithm = _default_best_algorithm(len(r))

	return self._original_compute(r, algorithm)

hypervolume._original_compute = hypervolume.compute
hypervolume.compute = _hypervolume_compute


def _hypervolume_exclusive(self, p_idx, r, algorithm = None):
	"""
	Compute the exlusive contribution to the total hypervolume by the point at index p_idx, given a reference point and the provided hypervolume algorithm (optional).

	USAGE:
		hv.exclusive(0, [5.0]*2)
		hv.exclusive(0, [5.0]*2, algorithm = hv_algorithm.optimal2d())
		* p_idx - index of the point
		* r - reference point
		* algorithm - optional argument: hypervolume algorithm used for the computation.
	"""
	if not isinstance(r, list):
		raise TypeError("Reference point must be a list of lists of real numbers, e.g.: r = [1.0, 1.0, 1.0]")
	try:
		r = [float(ri) for ri in r]
	except ValueError:
		raise TypeError("Every item in reference point (r), must be castable to float, e.g.: r = [1, '2.5', 10e-4]")

	if (algorithm == None):
		algorithm = _default_best_algorithm(len(r))

	return self._original_exclusive(p_idx, r, algorithm)

hypervolume._original_exclusive = hypervolume.exclusive
hypervolume.exclusive = _hypervolume_exclusive


def _hypervolume_least_contributor(self, r, algorithm = None):
	"""
	Find the least contributing point among the pareto front approximation.

	USAGE:
		hv.least_contributor([5.0]*2)
		hv.least_contributor([5.0]*2, algorithm = hv_algorithm.optimal2d())
		* r - reference point
		* algorithm - optional argument: hypervolume algorithm used for the computation.
	"""
	if not isinstance(r, list):
		raise TypeError("Reference point must be a list of lists of real numbers, e.g.: r = [1.0, 1.0, 1.0]")
	try:
		r = [float(ri) for ri in r]
	except ValueError:
		raise TypeError("Every item in reference point (r), must be castable to float, e.g.: r = [1, '2.5', 10e-4]")

	if (algorithm == None):
		algorithm = _default_best_algorithm(len(r))

	return self._original_least_contributor(r, algorithm)

hypervolume._original_least_contributor = hypervolume.least_contributor
hypervolume.least_contributor = _hypervolume_least_contributor

def _hypervolume_get_nadir_point(self, eps = 0.0):
	"""
	Return Nadir point for given set of points.

	USAGE:
		hv.nadir_point(eps = 1.0)
		* eps - value added to every objective in order to assert a strong dominance of reference point (0.0 by default).
	"""
	if eps < 0.0:
		raise TypeError("Epsilon must be a positive value.")

	return self._original_get_nadir_point(eps)

hypervolume._original_get_nadir_point = hypervolume.get_nadir_point
hypervolume.get_nadir_point = _hypervolume_get_nadir_point

def _optimal2d_ctor(self):
	"""
	Hypervolume algorithm: Optimal2D algorithm.
	Points are initially sorted by one dimension, after which the partial areas are summed linearly.
	Computational complexity: O(n*logn)

	USAGE:
		hv.compute(r=[1.5]*2, algorithm = hv_algorithm.optimal2d())
	"""
	return self._original_init()

optimal2d._original_init = optimal2d.__init__
optimal2d.__init__ = _optimal2d_ctor

def _optimal3d_ctor(self):
	"""
	Hypervolume algorithm: Optimal3D algorithm.
	Computational complexity: O(n*logn)

	REF: "On the Complexity of Computing the Hypervolume Indicator", Nicola Beume, Carlos M. Fonseca, Manuel Lopez-Ibanez,
	Luis Paquete, Jan Vahrenhold.
	IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTTATION. VOL. 13, NO. 5, OCTOBER 2009

	USAGE:
		hv.compute(r=[1.5]*3, algorithm = hv_algorithm.optimal3d())
	"""
	return self._original_init()

optimal3d._original_init = optimal3d.__init__
optimal3d.__init__ = _optimal3d_ctor

def _lebmeasure_ctor(self):
	"""
	Hypervolume algorithm: LebMeasure algorithm.
	Computational complexity: O(n*logn)

	REF: "A new analysis of the LebMeasure Algorithm for Calculating Hypervolume", L. While.

	USAGE:
		hv.compute(r=ref_point, algorithm = hv_algorithm.lebmeasure())
	"""
	return self._original_init()

lebmeasure._original_init = lebmeasure.__init__
lebmeasure.__init__ = _lebmeasure_ctor
