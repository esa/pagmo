from _base_stochastic import base_stochastic

class py_example_stochastic(base_stochastic):
	"""
	Noisy De Jong function implemented purely in Python.
	"""
	def __init__(self, dim = 1, seed = 0):
		super(py_example_stochastic,self).__init__(dim, seed)
		self.set_bounds(-5.12,5.12)
		self._dim = dim
	def _objfun_impl(self,x):
		from numpy.random import seed, rand	
		#We initialize the random number geneator of numpy using the 
		#data member seed. This will be changed by the algorithm whenever it detects 
		#a stochastic problem.
		seed(self.seed)
		print self.seed
		
		#And now we write the objfun that will always use the same pseudorandonm sequence
		#as long as self.seed is unchanged.
		f = 0;
		for i in range(self._dim):
			noise = rand()/10
			f = f + (x[i] + noise)*(x[i] + noise)
		return (f,)
