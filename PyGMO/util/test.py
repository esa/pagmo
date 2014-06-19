from PyGMO import *
from _analysis import analysis
from my_module import my_problem
from numpy import *
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

prob=problem.himmelblau()
anal=analysis(prob,1000, output_to_file=False)
# anal.f_distribution(percentile=[5,10,25,50,75],plot1=True,plot2=True,round_to=3)
# anal.f_linearity_convexity(n_pairs=0,tol=10**(-10),round_to=3)
# anal.f_regression(degree=[1,1,2,3],interaction=[True,False,False,False],pred=False,tol=10**(-8),round_to=3)
# anal.f_correlation(round_to=3)
# anal.c_linearity(npairs=0,tol=10**(-10),round_to=3)
# anal.c_feasibility(tol=10**(-8),round_to=3)
#anal.c_regression(degree=[1,2],interaction=False,pred=True,tol=10**(-8),round_to=3)
anal.multimodality(cluster=True,clusters_to_show=10, sample_size=0,algo=algorithm.gsl_fr(),decomposition_method='tchebycheff',\
    weights='uniform',z=[],variance_ratio=0.95,k=0,single_cluster_tolerance=0.001,kmax=0,round_to=3)

# anal2.f_distribution(percentile=[5,10,25,50,75],plot1=True,plot2=True,round_to=3)
# anal2.f_linearity_convexity()
# anal2.c_linearity()
# anal2.c_feasibility()
# anal.f_regression(degree=[1,2],interaction=[True,False],pred=True,tol=10**(-8),round_to=3)
# anal2.c_regression(degree=[1,1,2,3],interaction=[True,False,False,False],pred=False,tol=10**(-8),round_to=3)
# anal2.f_correlation()
# anal2.c_linearity()
# anal._get_gradient(mode='c')
# anal.plot_gradient_pcp(mode='c',invert=False)
# anal.plot_gradient_pcp(mode='c',invert=True)



