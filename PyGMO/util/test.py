from PyGMO import *
from _analysis import analysis
from my_module import my_problem
from numpy import *
import scipy as sp
import numpy as np

prob=problem.kur(2)
anal=analysis(prob)
anal.print_report()
