from PyGMO import *
from _analysis import analysis
from my_module import my_problem
from numpy import *
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

# prob=problem.himmelblau()
# anal=analysis(prob)
# anal.sample(1000)
# anal.get_local_extrema()
# anal.cluster_local_extrema()
# anal.plot_local_cluster_scatter()

prob=problem.kur(5)
anal=analysis(prob)
anal.sample(100)
anal.get_local_extrema(decomposition_method='bi',weights=[0.3,0.7],z=[0.1,0.1])
# anal.plot_local_cluster_scatter()
#anal.print_report(1000)