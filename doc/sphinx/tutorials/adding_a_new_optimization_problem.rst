.. _adding_a_new_optimization_problem:

================================================================
Adding a new optimization problem
================================================================

In this Tutorial we will learn how to code simple optimization problems
(continuous, single objective, unconstrained), so that PyGMO can then apply all of its
algorithmic power to solve it.

In a nutshell .... we will write classes deriving from problem.base
and reimplement some of its 'virtual' methods.

.. code-block:: python

   from PyGMO.problem import base
   class my_problem(base):
   """
   De Jong (sphere) function implemented purely in Python.
   
   USAGE: my_problem(dim = 10)

   * dim problem dimension
   """
   def __init__(self, dim = 10):
      #First we call the constructor of the base class telling
      #essentially to PyGMO what kind of problem to expect (1 objective, 0 contraints etc.)
      super(my_problem,self).__init__(dim)

      #then we set the problem bounds (in this case equal for all components)
      self.set_bounds(-5.12,5.12)
 
      #we define some additional 'private' data members (not really necessary in
      #this case, but ... hey this is a tutorial)
      self.__dim = dim

   #We reimplement the virtual method that defines the objective function.
   def _objfun_impl(self,x):
      f = 0;
      for i in range(self.__dim):
         f = f + (x[i])*(x[i])
         #note that we return a tuple with one element only. In PyGMO the objective functions
         #return tuples so that multi-objective optimization is also possible.
      return (f,)

   #Finally we also reimplement a virtual method that adds some output to the __repr__ method
   def human_readable_extra(self):
      return "\n\t Problem dimension: " + str(self.__dim)

Note that by default PyGMO will assume one wants to minimize the objective function. In the second
part of this tutorial we will also see how it is possible to change this default behaviour.

We may then put the above code in a file, say my_module.py and use, for example, Artificial Bee Colony .... with
20 individuals ....

.. code-block:: python

   from PyGMO import *
   from  my_module import my_problem

   prob = my_problem(dim=10)
   algo = algorithm.bee_colony(gen=500)
   isl = island(algo,prob,20)
   isl.evolve(1); isl.join()
   print isl.population.champion.f

And we are done!!!!! (the output will be something like 10^-27, no big deal for a sphere problem)

Let's consider now a maximization problem. To solve such a problem, two possibilities
are available to the PaGMO/PyGMO user. The first one is to code the original
problem as a minimization problem by premultiplying the objective function by -1 (a technique
wich is often used and requires no particular effort). If such a method is used,
the final fitness value obtained with PyGMO has to be multiplied by -1 to
get back to the correct value.

A second method, more elegant and most of all serving the purpose to show the use
of another virtual method which can be reimplemented in python objects deriving from base,
is to override the function that compares two fitness vectors. This function is used
by all pagmo algorithms to compare performances of individuals. By default, this function
compares the fitness f1 to a fitness f2 and returns true if f1 dominates f2 (which is single
objective optimization correspond to minimization). Let us see how ....

.. code-block:: python

    from PyGMO.problem import base
    class my_problem(base):
        """
        Analytical function to maximize
    
        USAGE: my_problem()
        """
        def __init__(self):
            # first we call the constructor of the base telling
            # to PyGMO what kind of problem to expect (1 objective, 0 constraints etc...)
            super(my_problem,self).__init__(2);
        
            # sets the problem bounds
            self.set_bounds(-10,10);
        
			# we do not need private members in this simple case
        
			# initialize best known solutions (this is optional and is here only
			# for demonstration purposes)
            self.best_x = [[1.,-1.]];
        
        # reimplement the virtual method that defines the obf function
        def _objfun_impl(self,x):
            f = ( - (1. - x[0])**2 - 100 * (-x[0]**2 - x[1])**2 - 1.);
            return(f,)
        
        # reimplement the virtual method that compares fitnesses
        def _compare_fitness_impl(self,f1,f2):
            return f1[0]>f2[0];

        # add some output to __repr__
        def human_readable_extra(self):
			return "\n\tMaximization problem"

As before, we may put this in the file my_module.py and use our favorite optimization
algorithm:

.. code-block:: python

    from PyGMO import *
    from my_module import my_problem
    from math import *

    prob = my_problem();
    algo = algorithm.de(gen=20);
    isl = island(algo,prob,20);
    isl.evolve(10); isl.join();

    print "Best individual:"
    print isl.population.champion

    print "Comparison of the best found fitness with the best known fitness:"
    for best_fitness in prob.best_f:
        print best_fitness[0] - isl.population.champion.f[0]

    print "L2 distance to the best decision vector"
    for best_decision in prob.best_x:
        l2_norm = 0;
        for n in range(0, len(best_decision)):
            l2_norm +=  (best_decision[n] - isl.population.champion.x[n])**2;
        l2_norm = sqrt(l2_norm);
        print l2_norm;

Note here that we used the best_f and best_x methods which return the best known
fitness and decision vectors. The best_f vector is automatically available as
we defined best_x in the problem. With these vectors, we can have an idea of
the optimizer performances. The result of this optimization is something
like 10^-11 for the comparison with the best fitness and 10^-5 for
the distance to the best decision vector.

NOTE1: This simple tutorial is implemented in PyGMO under the name PyGMO.problem.py_example
and PyGMO.problem.py_example_max

NOTE2: When evolve is called from an island, the process is forked and transferred to another python or ipython
instance. As a consequence, when writing your _obj_fun_impl you cannot use stuff like matplotlib to 
make interactive plots and alike. If you need, during development, to have this kind of support,
use the algorithm evolve method, for example

.. code-block:: python

    from PyGMO import *
    from  my_module import my_problem

    prob = my_problem(dim=10)
    algo = algorithm.bee_colony(gen=100)
    isl = island(algo,prob,20)
    pop = island.population
    pop = algo.evolve(pop)
    
    print "Best fitness:"
    print pop.champion.f

    print "Fitness found compared to the best known fitness:"
    for best_fitness in prob.get_best_known_f_vector():
        print pop.champion.f[0] - best_fitness[0]

NOTE3: If performance is your goal, you should implement your problem in C++, and then expose it into python.
