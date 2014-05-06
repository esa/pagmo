from PyGMO.problem import base
class my_problem(base):

   """De Jong (sphere) function implemented purely in Python.

   USAGE: my_problem(dim = 10)

   dim problem dimension
   """

   def __init__(self, dim = 10):
      #First we call the constructor of the base class telling
      #essentially to PyGMO what kind of problem to expect (1 objective, 0 contraints etc.)
      super(my_problem,self).__init__(dim)

      #then we set the problem bounds (in this case equal for all components)
      self.set_bounds(-5,5)

      #we define some additional 'private' data members (not really necessary in
      #this case, but ... hey this is a tutorial)
      self.__dim = dim
      #self.best_x=[[0,0,0,0,0,0,0,0,0,0]]

   #We reimplement the virtual method that defines the objective function.
   def _objfun_impl(self,x):
      f=20
      for i in range(self.__dim):
         f = f + i*(x[i])**2
         #note that we return a tuple with one element only. In PyGMO the objective functions
         #return tuples so that multi-objective optimization is also possible.
      s=(f,)
      return s

   #Finally we also reimplement a virtual method that adds some output to the __repr__ method
   def human_readable_extra(self):
      return "\n\t Problem dimension: " + str(self.__dim)