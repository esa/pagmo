from PyGMO import *
from _analysis import analysis
from my_module import my_problem
from numpy import *
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

prob=problem.cec2006(1)
anal=analysis(prob,1000,output_to_file=False,first=0)
anal._compute_constraints()

print anal._c_lin()