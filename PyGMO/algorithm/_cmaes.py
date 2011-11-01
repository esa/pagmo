from _base import base

class py_cmaes(base):
	"""
	Covariance Matrix Adaptation Evolutionary Strategy (Python)
	"""
	def __init__(self, gen = 500, elite = 0.5, cc = -1, cs = -1, c1 = -1, cmu = -1, sigma0=0.5, screen_output = False):
		"""
		Covariance Matrix Adaptation Evolutionary Strategy (Python)

		USAGE: algorithm.py_cmaes(gen = 500, elite = 0.5, cc = -1, cs = -1, c1 = -1, cmu = -1, screen_output = False)

		NOTE: In our variant of the algorithm, particle memory is used to extract the elite and reinsertion
		is made aggressively ..... geting rid ofthe worst guy). Also the bounds of the problem
		are enforced as to allow PaGMO machinery to work

		* gen: number of generations
		* elite: fraction of the population considered as elite (in (0,1])
		* cc: time constant for C cumulation (in [0,1]) if -1 automatic values are set
		* cs: time constant for sigma cumulation (in [0,1]) if -1 automatic values are set
		* c1: learning rate for rank-1 update (in [0,1]) if -1 automatic values are set
		* cmu: learning rate for rank-mu update (in [0,1]) if -1 automatic values are set
		* sigma0: starting step (std)
		* screen_output: activates screen_output (output at each generation)
		"""
		try:
			import numpy as np
		except ImportError:
			raise ImportError("This algorithm needs numpy to run. Is numpy installed?")

		if ( (cc<0 or cc>1) and not cc==-1):
			raise ValueError("cc needs to be in [0,1]")

		if ( (cs<0 or cs>1) and not cc==-1):
			raise ValueError("cs needs to be in [0,1]")

		if ( (c1<0 or c1>1) and not cc==-1):
			raise ValueError("c1 needs to be in [0,1]")

		if ( (cmu<0 or cmu>1) and not cc==-1):
			raise ValueError("cmu needs to be in [0,1]")

		base.__init__(self)
		self.__cc = cc
		self.__cs = cs
		self.__c1 = c1
		self.__cmu = cmu
		self.__gen = gen
		self.__elite = elite
		self.__sigma0 = sigma0
		self.screen_output = screen_output
		np.random.seed()

	def evolve(self,pop):
		from numpy import matrix,array,log, diag, eye,sqrt, exp, ones
		from numpy.random import normal, random
		from numpy.linalg import norm, eig

		# Let's rename some variables
		prob = pop.problem
		lb = prob.lb
		ub = prob.ub
		dim, cont_dim, int_dim, c_dim = prob.dimension, prob.dimension - prob.i_dimension, prob.i_dimension, prob.c_dimension

		# And perform checks on the problem type
		if cont_dim == 0:
			raise ValueError("There is no continuous dimension for cross_entropy to optimise!!")

		if c_dim > 0:
			raise ValueError("This version of cross_entropy is not suitable for constrained optimisation")

		if int_dim > 0:
			raise ValueError("The chromosome has an integer part .... this version of cross_entropy is not able to deal with it")

		# Setting sizes .....
		N = dim
		lam = len(pop)
		mu = lam/2

		# Setting coefficients for Selection
		weights = [log(mu+0.5) - log(i+1) for i in range(mu)]
		weights = [w / sum(weights) for w in weights]			# weights for weighted recombination
		mueff = sum(weights)**2 / sum(w**2 for w in weights)		# variance-effectiveness of sum w_i x_i

		# Setting coefficients for Adaptation
		if self.__cc == -1:
			cc = (4 + mueff/N) / (N+4 + 2*mueff/N);			# t-const for cumulation for C
		if self.__cs == -1:
			cs = (mueff+2) / (N+mueff+5);				# t-const for cumulation for sigma control
		if self.__c1 == -1:
			c1 = 2 / ((N+1.3)**2+mueff);				# learning rate for rank-one update of C
		if self.__cmu == -1:
			cmu = 2 * (mueff-2+1/mueff) / ((N+2)**2+mueff);		# and for rank-mu update

		damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs;		#damping for sigma

		# Initializing and allocationg
		mean = matrix(pop.champion.x).T
		variation = array([[0.0]*N]*lam)
		newpop = matrix([[0.0]*lam]*N)
		chiN = N**0.5*(1-1/(4*N)+1/(21*N**2))				#expectation of ||N(0,I)|| == norm(randn(N,1))
		B = matrix(eye(N,N));						#B defines the coordinate system
		D = ones(N);							#diagonal D defines the scaling
		C = matrix(eye(N,N));						#covariance matrix C
		invsqrtC = matrix(eye(N,N));					#inverse of sqrt(C)
		pc =matrix([[0]]*N)
		ps = matrix([[0]]*N)
		counteval = 0
		eigeneval = 0
		sigma=self.__sigma0

		# Let's start the algorithm
		for gen in range(self.__gen):

			#1 - We generate and evaluate lam new individuals
			variation = [B * diag(D) * normal(0,1,[dim,1]) for i in range(lam)]
			variation = [[j[0,0] for j in matr] for matr in variation]
			for i,d_mu in enumerate(variation):
				newpop[:,i] = mean + sigma * matrix(d_mu).T

			#fixing the bounds
			for row in range(newpop.shape[0]):
				for col in range(newpop.shape[1]):
					if newpop[row,col] > ub[row]:
						newpop[row,col] = lb[row] + random()*(ub[row]-lb[row])
					elif newpop[row,col] < lb[row]:
						newpop[row,col] = lb[row] + random()*(ub[row]-lb[row])

			#insert in population
			for i in range(lam):
				idx = pop.get_worst_idx()
				pop.set_x(idx,[newpop[j,i] for j in range(N)])
			counteval += lam

			#2 - We extract the elite from this generation
			#a = sorted(pop,lambda x,y: cmp(x.cur_f,y.cur_f))
			elite = [matrix(pop[idx].best_x).T for idx in pop.get_best_idx(mu)]
			#elite = [matrix(ind.cur_x).T for ind in a]
			#elite = elite[:mu]

			#3 - We compute the new elite mean storing the ol one
			meanold=mean
			mean = elite[0]*weights[0]
			for i in range(1,mu):
				mean += elite[i]*weights[i]

			#4 Update evolution paths
			ps = (1 - cs)*ps + sqrt(cs*(2-cs)*mueff)* invsqrtC * (mean-meanold) / sigma
			hsig = ((ps.T*ps)[0,0] / (1-(1-cs)**(2*counteval/lam)) / N) < (2 + 4/(N+1));
			hsig = int(hsig)
			pc = (1-cc) * pc + hsig * sqrt(cc*(2-cc)*mueff) * (mean-meanold) / sigma;

			#5 - Adapt Covariance Matrix
			Cold = C
			C = (elite[0]-meanold)*(elite[0]-meanold).T*weights[0]
			for i in range(1,mu):
				C += (elite[i]-meanold)*(elite[i]-meanold).T*weights[i]
			C /= sigma**2
			C = (1-c1-cmu)*Cold + cmu*C + c1 * ((pc * pc.T) + (1-hsig) * cc*(2-cc) * Cold)

			#6 - Adapt sigma
			sigma *= exp( (cs/damps)*(norm(ps)/chiN - 1));

			#7 - Perform eigen-decomposition of C
			#if ( (counteval - eigeneval) > (lam/(c1+cmu)/N/10) ):		#achieve O(N^2)
			eigeneval = counteval;
			C = (C+C.T)/2					#enforce symmetry
			D,B = eig(C);					#eigen decomposition, B==normalized eigenvectors
			D = [s**0.5 for s in D]				#D contains standard deviations now
			#if not (0 in D):				#Avoids numerical nans skipping evaluation of invsqrtC
			invsqrtC = B*diag([1/d for d in D])*B.T


			#8 - We print to screen if necessary
			if self.screen_output:
				if not(gen%20):
					print "\nGen.\tChampion\tHighest\t\tLowest\t\tVariation\t\tstep"
				print "%d\t%e\t%e\t%e\t%e\t%e" % (gen,pop.champion.f[0],
				    max([ind.cur_f[0] for ind in pop]),min([ind.cur_f[0] for ind in pop]),
				    norm(d_mu), sigma)
		return pop


	def get_name(self):
		return "CMAES (Python)"
	def human_readable_extra(self):
		return "gen=" + str(self.__gen) + " elite fraction=" + str(self.__elite) + " cc=" + str(self.__cc) + " cs=" + str(self.__cs) + " c1=" + str(self.__c1) + " cmu=" + str(self.__cmu) + " sigma0=" + str(self.__sigma0)
