# -*- coding: utf-8 -*-
from _util import hypervolume, hv_algorithm
from _util.hv_algorithm import lebmeasure, native2d, beume3d, wfg
from ..core._core import population

__all__ = ['hypervolume', 'hv_algorithm']

hv_algorithm.__doc__ = """Module containing available algorithms for the hypervolume computation

	USAGE:
		hv_algorithm.lebmeasure()
		hv_algorithm.native2d()
		hv_algorithm.beume3d
		hv_algorithm.wfg()
"""

class HypervolumeValidation:
	"""
	Utility class containing commonly raised errors.
	Kept in once place to simplify the consistency of error messages across methods
	"""

	# Raised when the reference point type is not a list or a tuple, e.g. r = "Foo"
	err_rp_type = TypeError("Reference point must be a list/tuple of real numbers, e.g.: r = [1.0, 1.0, 1.0]")

	# Raised when the reference point is a tuple/list but the items are non-castable to float, e.g. r = [1.0, 2.0, 'foo']
	err_rp_items_type = TypeError("Every item in reference point list/tuple must be castable to float, e.g.: r = [1, '2.5', 10e-4]")

	# Raised when the hypervolume object is constructed by anything other than a population object, tuple or a list, e.g. hypervolume("foo bar")
	err_hv_ctor_type = TypeError("Hypervolume object must be constructed from a list/tuple of points or a population object")

	@classmethod
	def handle_refpoint(cls, hypvol, r):
		"""
		Common way of validation for the reference point being passed as parameter to 'compute', 'exclusive' and 'least_contributor' methods.
		1. Check if user provided the reference point (use nadir point otherwise)
		2. Make sure that the reference point is of correct type
		3. Make sure that items of the reference point vector are castable to float
		"""
		if r:
			allowed_types = (list, tuple,)
			if not any(isinstance(r, T) for T in allowed_types):
				raise cls.err_rp_type
			try:
				r = [float(ri) for ri in r]
			except ValueError:
				raise cls.err_rp_items_type
		else:
			r = hypvol.get_nadir_point()
		return r

	@classmethod
	def handle_algorithm(cls, algorithm):
		"""
		Common way of validation for the hv_algorithm object being passed as parameter to 'compute', 'exclusive' and 'least_contributor' methods.
		"""
		hv_algo_types = (native2d, beume3d, wfg, lebmeasure, )
		if not any(isinstance(algorithm, T) for T in hv_algo_types):
			raise TypeError("'algorithm' keyword must be an instance of a hypervolume algorithm, e.g.: algorithm=hv_algorithm.wfg()")
		return algorithm


def _hypervolume_ctor(self, data_src = None, *args, **kwargs):
	"""
	Constructs a hypervolume object used for the computation of hypervolue and exclusive hypervolume.

	Object can be constructed from the population object, or from a fixed list/tuple of points
	Points within a fixed list must all be of equal size, of dimension larger than 1.

	USAGE:
		from PyGMO import *
		hv = hypervolume(pop)
		hv = hypervolume([[1,1,2],[2,1,2],[2,2,3]])
		hv = hypervolume(((1,2), (3,0.5), (1.5, 1.5))

	"""
	if not data_src or len(args) > 0 or len(kwargs) > 0:
		raise TypeError("Hypervolume takes exactly one argument in the constructor: population, list or a tuple object")

	allowed_types = (population, list, tuple,)
	if not any(isinstance(data_src, T) for T in allowed_types):
		raise HypervolumeValidation.err_hv_ctor_type

	args = []
	args.append(data_src)
	try:
		return self._original_init(*args)
	except TypeError:
		raise HypervolumeValidation.err_hv_ctor_type
hypervolume._original_init = hypervolume.__init__
hypervolume.__init__ = _hypervolume_ctor

def _hypervolume_compute(self, r = None, algorithm = None, *args, **kwargs):
	"""
	Compute the hypervolume indicator for a given reference point, using the provided hypervolume algorithm.
	Type 'hv_algorithm?' for a list of available hypervolume algorithms.

	USAGE:
		hv.compute()
		hv.compute(r=[5.0]*2)
		hv.compute(r=[5.0]*2, algorithm = hv_algorithm.native2d())
		* r (optional) - reference point, uses nadir point with default epsilon (see 'hypervolume.get_nadir_point?')
		* algorithm (optional) - hypervolume algorithm used for the computation, uses the best performing algorithm for given dimension by default.
	"""
	if len(args) > 0 or len(kwargs) > 0:
		raise TypeError("Incorrect combination of args/kwargs, type 'hypervolume.compute?' for usage")

	r = HypervolumeValidation.handle_refpoint(self, r)
	args = []
	args.append(r)
	if algorithm:
		algorithm = HypervolumeValidation.handle_algorithm(algorithm)
		args.append(algorithm)
	return self._original_compute(*args)

hypervolume._original_compute = hypervolume.compute
hypervolume.compute = _hypervolume_compute

def _hypervolume_exclusive(self, p_idx = None, r = None, algorithm = None, *args, **kwargs):
	"""
	Compute the exlusive contribution to the total hypervolume by the point at index p_idx, given a reference point and the provided hypervolume algorithm.
	Type 'hv_algorithm?' for a list of available hypervolume algorithms.

	USAGE:
		hv.exclusive(p_idx=0)
		hv.exclusive(p_idx=0, r=[5.0]*2)
		hv.exclusive(p_idx=0, r=[5.0]*2, algorithm=hv_algorithm.native2d())
		* p_idx - index of the point
		* r (optional) - reference point, uses nadir point with default epsilon (see 'hypervolume.get_nadir_point?')
		* algorithm (optional) - hypervolume algorithm used for the computation, uses the best performing algorithm for given dimension by default
	"""
	if p_idx == None:
		raise TypeError("p_idx (non-negative integer) argument is required for computation, type 'hypervolume.exclusive?' for usage.")
	if len(args) > 0 or len(kwargs) > 0:
		raise TypeError("Incorrect combination of args/kwargs, type 'hypervolume.exclusive?' for usage.")

	if not isinstance(p_idx, int) or p_idx < 0:
		raise TypeError("individual index (p_idx) must be a non-negative integer")

	r = HypervolumeValidation.handle_refpoint(self, r)

	args = []
	args.append(p_idx)
	args.append(r)
	if algorithm:
		algorithm = HypervolumeValidation.handle_algorithm(algorithm)
		args.append(algorithm)
	return self._original_exclusive(*args)

hypervolume._original_exclusive = hypervolume.exclusive
hypervolume.exclusive = _hypervolume_exclusive

def _hypervolume_least_contributor(self, r = None, algorithm = None, *args, **kwargs):
	"""
	Find the least contributing point among the pareto front approximation.
	Type 'hv_algorithm?' for a list of available hypervolume algorithms.

	USAGE:
		hv.least_contributor()
		hv.least_contributor(r=[5.0]*3)
		hv.least_contributor(r=[5.0]*3, algorithm=hv_algorithm.beume3d())
		* r (optional) - reference point, uses nadir point with default epsilon (see 'hypervolume.get_nadir_point?')
		* algorithm (optional) - hypervolume algorithm used for the computation, uses the best performing algorithm for given dimension by default
	"""

	if len(args) > 0 or len(kwargs) > 0:
		raise TypeError("Incorrect combination of args/kwargs, type 'hypervolume.least_contributor?' for usage")
	r = HypervolumeValidation.handle_refpoint(self, r)
	args = []
	args.append(r)
	if algorithm:
		algorithm = HypervolumeValidation.handle_algorithm(algorithm)
		args.append(algorithm)
	return self._original_least_contributor(*args)

hypervolume._original_least_contributor = hypervolume.least_contributor
hypervolume.least_contributor = _hypervolume_least_contributor

def _hypervolume_get_nadir_point(self, eps = 1.0):
	"""
	Return Nadir point for given set of points.

	USAGE:
		hv.nadir_point(eps = 1.0)
		* eps (optional) - value added to every objective in order to assert a strong dominance of reference point (1.0 by default).
	"""
	if eps < 0.0:
		raise ValueError("Epsilon must be a positive value.")

	return self._original_get_nadir_point(eps)

hypervolume._original_get_nadir_point = hypervolume.get_nadir_point
hypervolume.get_nadir_point = _hypervolume_get_nadir_point

def _native2d_ctor(self):
	"""
	Hypervolume algorithm: Native2D.
	Points are initially sorted by one dimension, after which the partial areas are summed linearly.
	Computational complexity: O(n*logn)
	Applicable to hypervolume computation problems of dimension=2

	USAGE:
		hv = hypervolume(...) # see 'hypervolume?' for usage
		refpoint=[1.0]*2
		hv.compute(r=refpoint, algorithm=hv_algorithm.native2d())
		hv.exclusive(p_idx=13, refpoint, algorithm=hv_algorithm.native2d())
		hv.least_contributor(r=refpoint, algorithm=hv_algorithm.native2d())
	"""
	return self._original_init()
native2d._original_init = native2d.__init__
native2d.__init__ = _native2d_ctor

def _beume3d_ctor(self):
	"""
	Hypervolume algorithm: Beume3D.
	Computational complexity: O(n*logn)
	Applicable to hypervolume computation problems of dimension=3

	REF: "On the Complexity of Computing the Hypervolume Indicator", Nicola Beume, Carlos M. Fonseca, Manuel Lopez-Ibanez,
	Luis Paquete, Jan Vahrenhold.
	IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTTATION. VOL. 13, NO. 5, OCTOBER 2009

	USAGE:
		hv = hypervolume(...) # see 'hypervolume?' for usage
		refpoint = [1.0]*3
		hv.compute(r=refpoint, algorithm=hv_algorithm.beume3d())
		hv.exclusive(p_idx=13, r=refpoint, algorithm=hv_algorithm.beume3d())
		hv.least_contributor(r=refpoint, algorithm=hv_algorithm.beume3d())
	"""
	return self._original_init()
beume3d._original_init = beume3d.__init__
beume3d.__init__ = _beume3d_ctor

def _lebmeasure_ctor(self):
	"""
	Hypervolume algorithm: LebMeasure.
	Computational complexity: O(n * 2^d) (very slow)
	Applicable to hypervolume computation problems of dimension in [3, ..]

	REF: "A new analysis of the LebMeasure Algorithm for Calculating Hypervolume", L. While.

	USAGE:
		hv = hypervolume(...) # see 'hypervolume?' for usage
		refpoint = [1.0]*7
		hv.compute(r=refpoint, algorithm=hv_algorithm.lebmeasure())
		hv.exclusive(p_idx=13, r=refpoint, algorithm=hv_algorithm.lebmeasure())
		hv.least_contributor(r=refpoint, algorithm=hv_algorithm.lebmeasure())
	"""
	return self._original_init()
lebmeasure._original_init = lebmeasure.__init__
lebmeasure.__init__ = _lebmeasure_ctor

def _wfg_ctor(self):
	"""
	Hypervolume algorithm: WFG.
	Applicable to hypervolume computation problems of dimension in [2, ..]

	REF: "A Fast Way of Calculating Exact Hypervolumes", Lyndon While, Lucas Bradstreet, Luigi Barone.
	IEEE TRANSACXTIONS ON EVOLUTIONARY COMPUTATION, VOL. 16, NO. 1, FEBRURARY 2012

	USAGE:
		hv = hypervolume(...) # see 'hypervolume?' for usage
		refpoint = [1.0]*7
		hv.compute(r=refpoint, algorithm=hv_algorithm.wfg())
		hv.exclusive(p_idx=13,r=refpoint, algorithm=hv_algorithm.wfg())
		hv.least_contributor(r=refpoint, algorithm=hv_algorithm.wfg())
	"""
	return self._original_init()
wfg._original_init = wfg.__init__
wfg.__init__ = _wfg_ctor
