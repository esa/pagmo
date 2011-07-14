# -*- coding: iso-8859-1 -*-
from _problem import *

# Raw C++ base class.
_base = _problem._base

# Import base python problem class
from _base_class import base

# Import PyGMO classes and methods imlemented in pure python (i.e. not exposed via boost)
from _spheres import visualize as _visualize
from _spheres_q import visualize as _visualize_q
spheres.visualize = _visualize
spheres_q.visualize = _visualize_q


from _pl2pl import py_pl2pl

class py_test(base):
	"""
	Minimal test problem implemented purely in Python. Objective function
	is De Jong's unidimensional sphere function.
	"""
	def __init__(self):
		super(py_test,self).__init__(1)
	def _objfun_impl(self,x):
		return (x[0] * x[0],)

def _get_problem_list():
	from PyGMO import problem
	return [problem.__dict__[n] for n in filter(lambda n: not n.startswith('_') and not n == 'base' and issubclass(problem.__dict__[n],problem._base),dir(problem))]





