# -*- coding: utf-8 -*-
from _util import hypervolume, hv_algorithm
from _util.hv_algorithm import lebmeasure, optimal2d, optimal3d

__all__ = ['hypervolume', 'hv_algorithm']

def _hypervolume_ctor(self, *args, **kwargs):
	"""
	Constructs a hypervolume object used for the computation of hypervolue and exclusive hypervolume.
	Constructors use unnamed arguments.

	USAGE:
		Construction by population object:

		hv = hypervolume(pop)
		* pop - instance of the population object

		Construction by fixed list of points:

		hv = hypervolume(points_list)

		* points_list - non-empty list of lists, all of equal size >= 2.
	"""
	return self._original_init(*args, **kwargs)

hypervolume._original_init = hypervolume.__init__
hypervolume.__init__ = _hypervolume_ctor

def _hypervolume_compute(self, r, algorithm = None):
	"""
	Compute the hypervolume indicator for a given reference point, using the provided hypervolume algorithm.

	USAGE:
		hv.compute(r=[5.0]*2)
		hv.compute(r=[5.0]*2, algorithm = hv_algorithm.optimal2d())
		* r - reference point
		* algorithm - optional argument: hypervolume algorithm used for computation. This method defaults to the best performing 
		algorithm for given dimension size.
	"""
	if not isinstance(r, list):
		raise TypeError("Reference point must be a list of lists of real numbers, e.g.: r = [1.0, 1.0, 1.0]")
	try:
		r = [float(ri) for ri in r]
	except ValueError:
		raise TypeError("Every item in reference point (r), must be castable to float, e.g.: r = [1, '2.5', 10e-4]")

	if (algorithm == None):
		if (len(r) == 2):
			algorithm = optimal2d()
		elif (len(r) == 3):
			algorithm = optimal3d()
		else:
			algorithm = lebmeasure()
	return self._original_compute(r, algorithm)

hypervolume._original_compute = hypervolume.compute
hypervolume.compute = _hypervolume_compute

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
