from PyGMO.problem import base


class py_example(base):

    """
    De Jong (sphere) function implemented purely in Python.

    USAGE: py_example(dim = 10)

    * dim problem dimension
    """

    def __init__(self, dim=10):
        # First we call the constructor of the base class telling
        # essentially to PyGMO what kind of problem to expect (1 objective, 0
        # contraints etc.)
        super(py_example, self).__init__(dim)

        # then we set the problem bounds (in this case equal for all
        # components)
        self.set_bounds(-5.12, 5.12)

        # and we define some additional 'private' data members (not really necessary in
        # this case, but ... hey this is a tutorial)
        self.__dim = dim

    # We reimplement the virtual method that defines the objective function.
    def _objfun_impl(self, x):
        f = 0
        for i in range(self.__dim):
            f = f + (x[i]) * (x[i])
        # note that we return a tuple with one element only. In PyGMO the objective functions
        # return tuples so that multi-objective optimization is also possible.
        return (f,)

    # Finally we also reimplement a virtual method that adds some output to
    # the __repr__ method
    def human_readable_extra(self):
        return "\n\t Problem dimension: " + str(self.__dim)


class py_example_max(base):

    """
    Analytical function to maximize

    USAGE: py_example_max()
    """

    def __init__(self):
        # first we call the constructor of the base telling
        # to PyGMO what kind of problem to expect (1 objective, 0 constraints
        # etc...)
        super(py_example_max, self).__init__(2)

        # sets the problem bounds
        self.set_bounds(-10, 10)

        # define private data members
        self.__dim = 2

        # initialize best known solutions
        self.best_x = [[1., -1.]]

    # reimplement the virtual method that defines the obf function
    def _objfun_impl(self, x):
        f = (- (1. - x[0]) ** 2 - 100 * (-x[0] ** 2 - x[1]) ** 2 - 1.)
        return(f,)

    # reimplement the virtual method that compares fitnesses
    def _compare_fitness_impl(self, f1, f2):
        return f1[0] > f2[0]

    # add some output to __repr__
    def human_readable_extra(self):
        return "\n\t Maximization problem"
