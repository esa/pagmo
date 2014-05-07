from PyGMO import *
from _analysis import analysis
from my_module import my_problem
from numpy import *
import scipy as sp

anal=analysis()
prob=problem.rastrigin(5)
anal.sample(prob,10000)
# # print anal.lda()
# # print anal.qda()
# # print anal.knn()
anal.get_local_extrema(prob)
print anal.local_extrema," \n"
print anal.local_f," \n\n"
# print anal.cluster_local_extrema()