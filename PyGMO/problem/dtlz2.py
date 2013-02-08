from PyGMO.problem import base
from mpl_toolkits.mplot3d import axes3d

import itertools
import matplotlib.pyplot as plt
import numpy as np

class dtlz2(base):
	"""
	DTLZ2 - scaling function test implementation.
	
	This problem is an easy problem for multi-objective optimization which can have an
	arbitraty large fitness dimension. The search space is continuos, unimodal and the
	problem is not deceptive.
	
	USAGE: prob = dtlz2(dim=10, fdim=3)
	
	* k dimension of the problem is k + fdim - 1
	* fdim fitness dimension
	"""
	def __init__(self, k = 10, fdim = 3):
		super(dtlz2,self).__init__(k + fdim - 1, 0, fdim)
		self.set_bounds(0.0,1.0)
		self.fdim = fdim
		self.k = k

	def dtlz2_g(self, x):
		""" The "g"-function of the dtlz2 problem."""
		return sum([(y - 0.5)**2 for y in x])

	def _get_typename(self):
		return 'DTLZ'
		
	def _objfun_impl(self,x):
		"""
		The chromosome: x_1, x_2, ..., ..., x_M-1, x_M, x_M+1, ... , x_M+k
												   [---( Vector x_M )-----]
					   x[0],x_[1], ..., x[fdim-2], x[fdim-1], ..., x[fdim-1+k]
		"""
		f = self.fdim * [0.0];
		g = self.dtlz2_g(x[self.fdim-1:])
		
		f[0] = (1 + g) * np.sin(x[0] * (np.pi/2))
		for idx, y in enumerate(f[1:-1]):
			if x[idx] == 0:
				f[idx+1] = 0			# otherwise we divide by zero
			else:
				f[idx+1] = f[idx] * np.cos(x[idx] * (np.pi/2)) * np.sin(x[idx+1] * (np.pi/2)) / np.sin(x[idx] * (np.pi/2))
		
		f[-1] = (1 + g) * np.product([np.cos(y * (np.pi/2)) for y in x[:self.fdim-1]])

		f.reverse()	# to be conform with the definition
		return tuple(f)
		
	def plot(self, pop, a=40):
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		
		# plot the wireframe of the known optimal pareto front
		thetas = np.linspace(0, (np.pi / 2.0), 30)
		#gammas = np.linspace(-np.pi / 4, np.pi / 4, 30)
		gammas = np.linspace(0, (np.pi / 2.0), 30)
		
		x_frame = np.outer(np.cos(thetas), np.cos(gammas))
		y_frame = np.outer(np.cos(thetas), np.sin(gammas))
		z_frame = np.outer(np.sin(thetas), np.ones(np.size(gammas)))
		
		ax.view_init(azim=a)
		
		ax.set_autoscalex_on(False)
		ax.set_autoscaley_on(False)
		ax.set_autoscalez_on(False)
		
		ax.set_xlim(0, 1.8)
		ax.set_ylim(0, 1.8)
		ax.set_zlim(0, 1.8)
		
		ax.plot_wireframe(x_frame,y_frame,z_frame)
		
		
		# plot the individuals of the population
		fit = np.transpose([ind.cur_f for ind in pop])
		ax.plot(fit[0],fit[1],fit[2], 'ro')
				
		
	def human_readable_extra(self):
		return "\n\t Problem dimension: " + str(self.k + self.fdim - 1)