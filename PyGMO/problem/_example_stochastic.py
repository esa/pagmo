from PyGMO.problem import base_stochastic


class py_example_stochastic(base_stochastic):

    """
    Noisy De Jong (sphere) function implemented purely in Python.

    USAGE: py_example_stochastic(dim = 10, seed=0)

    * dim problem dimension
    * seed initial random seed
    """

    def __init__(self, dim=10, seed=0):
        # First we call the constructor of the base stochastic class. (Only
        # unconstrained single objective problems can be stochastic in PyGMO)
        super(py_example_stochastic, self).__init__(dim, seed)

        # then we set the problem bounds (in this case equal for all
        # components)
        self.set_bounds(-5.12, 5.12)

        # and we define some additional 'private' data members (not really necessary in
        # this case, but ... hey this is a tutorial)
        self.__dim = dim

    def _objfun_impl(self, x):
        from random import random as drng
        from random import seed

        # We initialize the random number generator using the
        # data member seed (in base_stochastic). This will be changed by suitable
        # algorithms when a stochastic problem is used. The mod operation
        # avoids overflows

        seed(self.seed)

        # We write the objfun using the same pseudorandonm sequence
        # as long as self.seed is unchanged.
        f = 0
        for i in range(self.__dim):
            noise = (2 * drng() - 1) / 10
            f = f + (x[i] + noise) * (x[i] + noise)
        return (f,)

    def human_readable_extra(self):
        return "\n\tSeed: " + str(self.seed)
