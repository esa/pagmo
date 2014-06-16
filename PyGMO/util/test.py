from PyGMO import *
from _analysis import analysis
from my_module import my_problem
from numpy import *
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

prob=problem.cec2006(5)
a=util.faure(5,1)
anal=analysis(prob,1000, output_to_file=False)
anal.f_distribution(percentile=[5,10,25,50,75],plot1=True,plot2=True,round_to=3)
anal.f_linearity_convexity()
anal.test_constraints()