# -*- coding: utf-8 -*-
from _hypervolume import hypervolume, lebmeasure, optimal2d, optimal3d

def _hypervolume_ctor(self, *args, **kwargs):
	"""
	Constructs a hypervolume object used for the computation of hypervolue and exclusive hypervolume.
	Constructors use unnamed arguments.

	USAGE:
		Construction by population object:

		hypervolume(pop)
		* pop - instance of the population object

		Construction by fixed list of points:

		hypervolume(points_list)

		* points_list - non-empty list of lists, all of equal size >= 2.
	"""
	return self._original_init(*args, **kwargs)
hypervolume._original_init = hypervolume.__init__
hypervolume.__init__ = _hypervolume_ctor

def _hypervolume_compute(self, r = None, algorithm = None):
	"""
	Compute the hypervolume indicator for a given reference point, using the provided hypervolume algorithm.

	USAGE:
		hv.compute(r=[5.0]*2, hv_algorithm=optimal2d())
	"""
	if not isinstance(r, list):
		raise TypeError("Reference point must be a list of lists of real numbers, e.g.: r = [1.0, 1.0, 1.0]")
	try:
		r = [float(ri) for ri in r]
	except ValueError:
		raise TypeError("Every item in reference point (r), must be castable to float, e.g.: r = [1, '2.5', 10e-4]")

	if (algorithm == None):
		if (len(r[0]) == 2):
			algorithm = optimal2d()
		else:
			algorithm = lebmeasure()
	return self._original_compute(r, algorithm)

hypervolume._original_compute = hypervolume.compute
hypervolume.compute = _hypervolume_compute

def _optimal2d_ctor(self):
	"""
	Hypervolume algorithm: Optimal2D algorithm.

	USAGE:
		optimal2d()
	"""
	return self._original_init()

optimal2d._original_init = optimal2d.__init__
optimal2d.__init__ = _optimal2d_ctor

def _optimal3d_ctor(self):
	"""
	Hypervolume algorithm: Optimal3D algorithm.

	USAGE:
		optimal3d()
	"""
	return self._original_init()

optimal3d._original_init = optimal3d.__init__
optimal3d.__init__ = _optimal3d_ctor

def _lebmeasure_ctor(self):
	"""
	Hypervolume algorithm: LebMeasure algorithm.

	USAGE:
		lebmeasure()
	"""
	return self._original_init()

lebmeasure._original_init = lebmeasure.__init__
lebmeasure.__init__ = _lebmeasure_ctor
