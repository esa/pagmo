from base import base

class py_cross_entropy(base):
       """
       Cross-Entropy algorithm (Python)
       """
       def __init__(self, gen = 1, elite = 0.5, scale = 0.05):
	      """
              Constructs a Cross-Entropy Algorithm (Python)
		
              USAGE: algorithm.py_cross_entropy(gen = 1, elite = 0.5, scale = 0.05)
		
              NOTE: A multivariate normal distribution is used. 
		    The first sample is centered around the population champion. 
		    Covariance matrix and mean is evaluated using ind.best_x

              * gen: number of generations
              * elite: fraction of the population considered as elite (in (0,1])
              * scale: scaling factor for the estimated covariance matrix
              """
              try:
			import numpy as np
              except ImportError:
			raise ImportError("This algorithm needs numpy to run. Is numpy installed?")
			
              base.__init__(self)
              self.__gen = gen
              self.__elite = elite
              self.__scale = scale
              np.random.seed()

       def evolve(self,pop):
	       from numpy import matrix
	       
               # Let's rename some variables
               prob = pop.problem
               dim, cont_dim, int_dim, c_dim = prob.dimension, prob.dimension - prob.i_dimension, prob.i_dimension, prob.c_dimension
               
               # And perform checks on the problem type
               if cont_dim == 0:
                       raise ValueError("There is no continuous dimension for cross_entropy to optimise!!")

               if c_dim > 0:
                       raise ValueError("This version of cross_entropy is not suitable for constrained optimisation")

               if int_dim > 0:
                       raise ValueError("The chromosome has an integer part .... this version of cross_entropy is not able to deal with it")
 
               # We then check that the elite is not empty
               n_ind__elite = int(len(pop) * self.__elite)
	       if n_ind__elite == 0:
                       raise ValueError("Elite contains no individuals ..... maybe increase the elite parameter?")

               # If the incoming population is empty ... do nothing
               if len(pop) == 0:
                       return population
                       
               # Let's start the algorithm                  
               mu = matrix(pop.champion.x)
               
               for _ in range(self.__gen):
                       y = self.__extract_elite(pop,n_ind__elite)		#y = array, [[chrom],fitness] * n_ind__elite
                       C = self.__estimate_covariance(y,mu) * self.__scale      #C = matrix, D x D
                       mu = self.__calculate_mean(y)                            #mu = matrix, D x 1
                       self.__new_generation(pop,mu,C,prob.lb,prob.ub)
               return pop

       def __extract_elite(self,pop,N):
	       from numpy import matrix, array
               # We transform the population into an easier to manipulate pythonic structure
               x = [[matrix(ind.best_x), ind.best_f[0]] for ind in pop]
                                               #x[i][0]: i-th chromosome (matrix)
                                               #x[i][1]: i-th fitness (scalar)
               x = array(x)
               x = sorted(x,key=lambda row: row[1])
               return x[:N]

       def __estimate_covariance(self,y,mu):
               from numpy import matrix
               D = y[0][0].size
               C = matrix([[0]*D]*D)
               for i in range(len(y)):
                       C = C + (y[i][0]-mu).T*(y[i][0]-mu)
               return (C / len(y))

       def __calculate_mean(self,x):
	       from numpy import matrix
               D = x[0][0].size
               mu = matrix([[0]*D])
               for i in range(len(x)):
                       mu = mu + x[i][0]
               return (mu.T / len(x))

       def __new_generation(self,pop,mu,C,lb,ub):
	       from numpy.random import multivariate_normal,random
	       from numpy import array
               np = len(pop)
               newpop = multivariate_normal(array(mu.T)[0],C,[np])
               for row in range(newpop.shape[0]):
                       for col in range(newpop.shape[1]):
                               if newpop[row,col] > ub[col]:
                                       newpop[row,col] = lb[col] + random()*(ub[col]-lb[col])
                               elif newpop[row,col] < lb[col]:
                                       newpop[row,col] = lb[col] + random()*(ub[col]-lb[col])
               for i in range(np):
                       pop.set_x(i,newpop[i])
               return
	

       def get_name(self):
               return "Cross Entropy (Python)"
       def human_readable_extra(self):
               return "gen=" + str(self.__gen) + " elite fraction=" + str(self.__elite) + " covariance scaling=" + str(self.__scale)
