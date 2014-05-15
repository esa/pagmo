from PyGMO import *
from _analysis import analysis
from my_module import my_problem
from numpy import *
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

prob=problem.kur(3)
anal=analysis(prob)
anal.print_report(1000)
