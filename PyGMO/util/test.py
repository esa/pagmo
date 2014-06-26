from PyGMO import *
from numpy import *
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

prob=problem.kur(4)

inspector=util.analysis(prob, 1000, output_to_file=True)
inspector.f_distribution()
# inspector.f_distribution(percentile=[10,25,50],plot_f_distribution=True,plot_x_pcp=True,round_to=3)

# inspector.f_linearity_convexity(n_pairs=0,tol=10**(-8),round_to=3)

# inspector.f_regression(degree=[1,1,2,2,3,3],interaction=[False,True,False,True,False,True],pred=[False,False,True,False,False,False],tol=10**(-8),round_to=3)

# inspector.f_sensitivity(hessian=True,plot_gradient_sparsity=True, plot_pcp=True, plot_inverted_pcp=True, \
#     sample_size=0,h=0.01,conv_tol=0.000001,zero_tol=0.000001,tmax=15,round_to=3)

# inspector.levelset(threshold=45,k_tune=2,k_test=10,linear=True,quadratic=True,nonlinear=True,round_to=3)

# inspector.local_search(clusters_to_show=10,plot_global_pcp=False,plot_separate_pcp=False,scatter_plot_dimensions=[],\
#     sample_size=0,algo=algorithm.gsl_fr(),decomposition_method='tchebycheff',weights='uniform',z=[],\
#     con2mo='obj_cstrsvio',variance_ratio=0.9,k=0,single_cluster_tolerance=0.001,kmax=0,round_to=3)

# inspector.f_correlation(tc=0.95,tabs=0.1,round_to=3)

# inspector.c_feasibility(tol=10**(-8),round_to=3)

# inspector.c_linearity(npairs=0,tol=10**(-8),round_to=3)

# inspector.c_regression(degree=[1,1,2,2],interaction=[False,True,False,True],pred=True,tol=10**(-8),round_to=3)

# inspector.c_sensitivity(plot_gradient_sparsity=True, plot_pcp=True, plot_inverted_pcp=True,sample_size=0,h=0.01,conv_tol=0.000001,zero_tol=0.000001,tmax=15,round_to=3)
