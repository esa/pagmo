from _base import base

class py_example(base):
	"""
	Shifted De Jong function implemented purely in Python.
	"""
	def __init__(self, dim = 1):
		super(py_example_stochastic,self).__init__(dim)
		self.set_bounds(-5.12,5.12)
		self._dim = dim
	def _objfun_impl(self,x):
		f = 0;
		for i in range(self._dim):
			f = f + (x[i])*(x[i])
		return (f,)
