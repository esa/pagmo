from PyGMO import *
from _analysis import analysis
from my_module import my_problem
from numpy import *
import scipy as sp

an=analysis()
prob=problem.kur(5)
an.sample(prob,1000,'lhs')
#an.compute_constraints(prob)
an.get_gradient(prob)
print an.grad_sparsity