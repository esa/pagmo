from PyGMO.problem._problem import _base


class base(_base):

    """
    This class is the base class for all non-stochastic optimization problems. When defining an optimization problem
    in PyGMO, the user needs to write a class that inherits from this base class and needs to call its constructor.
    He will then need to re-implement a number of virtual functions that define the problem objectives and constraints,
    as well as defining the box-bounds on the decision vector.
    """

    def __init__(self, *args):
        """
        Base problem constructor. It must be called from within the derived class constructor __init__()

        USAGE: super(derived_class_name,self).__init__(dim, i_dim, n_obj, c_dim, c_ineq_dim, c_tol)

        * dim: Total dimension of the decision vector
        * i_dim: dimension of the integer part of decision vector (the integer part is placed at the end of the decision vector). Defaults to 0
        * n_obj: number of objectives. Defaults to 1
        * c_dim: total dimension of the constraint vector. dDefaults to 0
        * c_ineq_dim: dimension of the inequality part of the constraint vector (inequality const. are placed at the end of the decision vector). Defaults to 0
        * c_tol: constraint tolerance. When comparing individuals, this tolerance is used to decide whether a constraint is considered satisfied.
        """
        if len(args) == 0:
            raise ValueError(
                "Cannot initialise base problem without parameters for the constructor.")
        _base.__init__(self, *args)

    def _get_typename(self):
        return str(type(self))

    def __get_deepcopy__(self):
        from copy import deepcopy
        return deepcopy(self)

    def get_name(self):
        return self._get_typename()
