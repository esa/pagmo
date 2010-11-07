# -*- coding: iso-8859-1 -*-
from _problem import *

class base(_problem._base):
	def __init__(self,*args):
		if len(args) == 0:
			raise ValueError("Cannot initialise base problem without parameters for the constructor.")
		_problem._base.__init__(self,*args)
	def _get_typename(self):
		return str(type(self))
	def get_name(self):
		return self._get_typename()

class py_test(base):
	def __init__(self):
		super(py_test,self).__init__(1)
	def __copy__(self):
		return py_test()
	def _objfun_impl(self,x):
		return (x[0] * x[0],)

class py_broken(base):
	def __init__(self):
		super(py_broken,self).__init__(1)
	def __copy__(self):
		return py_broken()
	def _objfun_impl(self,x):
		return (float('NaN'),)
