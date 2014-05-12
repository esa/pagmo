from PyGMO import *
from _analysis import analysis
from my_module import my_problem
from numpy import *
import scipy as sp
import numpy as np

anal=analysis()
prob=problem.schwefel(1)
anal.sample(prob,500)
anal.get_local_extrema(prob)
anal.cluster_local_extrema()
anal.print_report(prob)