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
