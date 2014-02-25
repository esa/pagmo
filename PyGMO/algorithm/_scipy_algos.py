from PyGMO.algorithm._base import base

# Helper class to ease the inclusion of scipy.optimize solvers.


class _scipy_base(base):

    def __init__(self, solver_name, constrained):
        base.__init__(self)
        try:
            exec('from scipy.optimize import %s as solver' % solver_name)
            from numpy import concatenate, array
        except ImportError:
            raise ImportError(
                'The necessary SciPy and/or NumPy classes/functions could not be imported')
        self.solver = solver
        self.constrained = constrained
    # Check if problem is compatible with the algorithm.

    def _problem_checks(self, prob):
        if prob.f_dimension > 1:
            raise ValueError(
                "this algorithm does not support multi-objective optimisation")
        if prob.dimension == prob.i_dimension:
            raise ValueError("the provided problem has no continuous part")
        if not self.constrained and prob.c_dimension:
            raise ValueError(
                "this algorithm does not support constrained optimisation")
    # Check that the algorithm did not go out of bounds, and, in such a case,
    # correct the chromosome.

    def _check_new_chromosome(self, new_chromosome, prob):
        for i in range(0, len(new_chromosome)):
            if new_chromosome[i] < prob.lb[i]:
                new_chromosome[i] = prob.lb[i]
            if new_chromosome[i] > prob.ub[i]:
                new_chromosome[i] = prob.ub[i]
        return new_chromosome

    def _starting_params(self, pop):
        from numpy import array
        # Number of equality constraints.
        n_ec = pop.problem.c_dimension - pop.problem.ic_dimension
        # Extract the continuous part of the first individual's current
        # chromosome.
        x0 = array(
            pop[pop.get_best_idx()].cur_x
            [0: pop.problem.dimension - pop.problem.i_dimension], dtype=float)
        # Combinatorial part of the chromosome (which will not be optimised).
        x0_comb = array(
            pop[pop.get_best_idx()].cur_x
            [pop.problem.dimension - pop.problem.i_dimension:], dtype=float)
        return n_ec, x0, x0_comb


class scipy_fmin(_scipy_base):

    """
    Wrapper around SciPy's fmin optimiser (Uses a Nelder-Mead simplex algorithm to find the minimum of function of one or more variables.)
    """

    def __init__(
            self,
            maxiter=1,
            xtol=0.0001,
            ftol=0.0001,
            maxfun=None,
            full_output=0,
            disp=0,
            retall=0):
        """
        Constructs a Nelder-Mead Simplex algorithm (SciPy)

        USAGE: algorithm.scipy_fmin(maxiter=1, xtol=0.0001, ftol=0.0001, maxfun=None, full_output=0, disp=0, retall=0)

        * maxiter: Maximum number of iterations to perform
        * xtol: Relative error in xopt acceptable for convergence
        * ftol: Relative error in func(xopt) acceptable for convergence
        * maxfun: Maximum number of function evaluations to make
        * full_output: Set to True if fval and warnflag outputs are desired
        * disp: Set to True to print convergence messages
        * retall: Set to True to return list of solutions at each iteration
        """
        _scipy_base.__init__(self, 'fmin', False)
        self.xtol = xtol
        self.ftol = ftol
        self.maxiter = maxiter
        self.maxfun = maxfun
        self.full_output = full_output
        self.disp = disp
        self.retall = retall

    def evolve(self, pop):
        from numpy import concatenate
        prob = pop.problem
        self._problem_checks(prob)
        if len(pop) == 0:
            return pop
        _, x0, x0_comb = self._starting_params(pop)
        retval = self.solver(
            lambda x: prob.objfun(
                concatenate(
                    (x,
                     x0_comb)))[0],
            x0,
            xtol=self.xtol,
            ftol=self.ftol,
            maxiter=self.maxiter,
            maxfun=self.maxfun,
            full_output=self.full_output,
            disp=self.disp,
            retall=self.retall)
        new_chromosome = list(retval) + list(x0_comb)
        pop.set_x(
            pop.get_best_idx(),
            self._check_new_chromosome(
                new_chromosome,
                prob))
        return pop

    def get_name(self):
        return "Nelder-Mead Simplex (SciPy)"

    def human_readable_extra(self):
        return "maxiter = " + str(self.maxiter) + ", xtol = " + str(self.xtol) + \
            ", ftol = " + str(self.ftol) + ", maxfun = " + str(self.maxfun)


class scipy_l_bfgs_b(_scipy_base):

    """
    Wrapper around SciPy's fmin_l_bfgs_b optimiser (uses L-BFGS-B algorithm)
    """

    def __init__(
            self,
            maxfun=1,
            m=10,
            factr=10000000.0,
            pgtol=1e-05,
            epsilon=1e-08,
            screen_output=False):
        """
        Constructs a L-BFGS-B algorithm (SciPy)

        NOTE: gradient is numerically approximated

        USAGE: algorithm.scipy_l_bfgs_b(maxfun = 15000, m = 10, factr = 10000000.0, pgtol = 1e-05, epsilon = 1e-08, screen_output = False):

        * maxfun: maximum number of function evaluations
        * m: the maximum number of variable metric corrections
                used to define the limited memory matrix. (the limited memory BFGS
                method does not store the full hessian but uses this many terms in an
                approximation to it).
        * factr: The iteration stops when
                (f{k} - f{k+1}) / max{\| f{k} \| , \| f{k+1} \|,1} <= factr*epsmch
                where epsmch is the machine precision, which is automatically
                generated by the code. Typical values for factr: 1e12 for
                low accuracy; 1e7 for moderate accuracy; 10.0 for extremely
                high accuracy.
        * pgtol: The iteration will stop when
                max{\| proj g{i} \| i = 1, ..., n} <= pgtol
                where proj g{i} is the ith component of the projected gradient.
        * epsilon: step size used when approx_grad is true, for numerically
                calculating the gradient
        * screen_output: Set to True to print iterations
        """
        _scipy_base.__init__(self, 'fmin_l_bfgs_b', False)
        self.maxfun = maxfun
        self.m = m
        self.factr = factr
        self.pgtol = pgtol
        self.epsilon = epsilon
        self.screen_output = screen_output

    def evolve(self, pop):
        from numpy import concatenate, array
        prob = pop.problem
        self._problem_checks(prob)
        if len(pop) == 0:
            return pop
        _, x0, x0_comb = self._starting_params(pop)
        # Extract the a list of tuples representing the bounds.
        prob_bounds = [(prob.lb[i], prob.ub[i])
                       for i in range(0, prob.dimension - prob.i_dimension)]
        if self.screen_output:
            iprn = 1
        else:
            iprn = -1
        retval = self.solver(
            lambda x: array(
                prob.objfun(
                    concatenate(
                        (x,
                         x0_comb))),
                dtype=float),
            x0,
            bounds=prob_bounds,
            approx_grad=True,
            iprint=iprn,
            pgtol=self.pgtol,
            maxfun=self.maxfun,
            factr=self.factr,
            m=self.m,
            epsilon=self.epsilon)
        new_chromosome = list(retval[0]) + list(x0_comb)
        pop.set_x(
            pop.get_best_idx(),
            self._check_new_chromosome(
                new_chromosome,
                prob))
        return pop

    def get_name(self):
        return "L-BFGS-B (SciPy)"

    def human_readable_extra(self):
        return "maxfun = " + str(self.maxfun) + ", m = " + str(self.m) + ", factr = " + \
            str(self.factr) + ", pgtol = " + str(self.pgtol) + ", epsilon = " + str(self.epsilon)


class scipy_slsqp(_scipy_base):

    """
    Wrapper around SciPy's slsqp optimiser.
    """

    def __init__(
            self,
            max_iter=100,
            acc=1E-8,
            epsilon=1.4901161193847656e-08,
            screen_output=False):
        """
        Constructs a Sequential Least SQuares Programming algorithm

        NOTE: gradient is numerically approximated

        USAGE: algorithm.scipy_slsqp(max_iter = 100,acc = 1E-6,epsilon = 1.49e-08, screen_output = False))


        * max_iter: The maximum number of iterations.
        * acc: Requested accuracy.
        * epsilon: The step size for finite-difference derivative estimates.
        * screen_output: Set to True to print iterations
        """

        _scipy_base.__init__(self, 'fmin_slsqp', True)
        self.max_iter = max_iter
        self.acc = acc
        self.epsilon = epsilon
        self.screen_output = screen_output

    def get_name(self):
        return 'Sequential Least SQuares Programming (SciPy)'

    def evolve(self, pop):
        from numpy import concatenate, array
        prob = pop.problem
        self._problem_checks(prob)
        # If population is empty, just return input population.
        if len(pop) == 0:
            return pop
        # Get starting params.
        n_ec, x0, x0_comb = self._starting_params(pop)
        # Extract the a list of tuples representing the bounds.
        prob_bounds = [(prob.lb[i], prob.ub[i])
                       for i in range(0, prob.dimension - prob.i_dimension)]
        if self.screen_output:
            iprn = 2
        else:
            iprn = 0
        # Run the optimisation.
        retval = self.solver(
            lambda x: prob.objfun(
                concatenate(
                    (x, x0_comb)))[0], x0, f_eqcons=lambda x: array(
                prob.compute_constraints(
                    concatenate(
                        (x, x0_comb)))[
                    0:n_ec], dtype=float), f_ieqcons=lambda x: array(
                prob.compute_constraints(
                    concatenate(
                        (x, x0_comb)))[
                    n_ec:], dtype=float) * -1, bounds=prob_bounds, iprint=iprn, iter=self.max_iter, acc=self.acc, epsilon=self.epsilon)
        # Set the individual's chromosome in the population and return. Conserve the integer part from the
        # original individual.
        new_chromosome = list(retval) + list(x0_comb)
        pop.set_x(
            pop.get_best_idx(),
            self._check_new_chromosome(
                new_chromosome,
                prob))
        return pop

    def human_readable_extra(self):
        return "maxiter = " + \
            str(self.max_iter) + ", acc = " + str(self.acc) + ", epsilon = " + str(self.epsilon)


class scipy_tnc(_scipy_base):

    """
    Wrapper around SciPy's tnc optimiser.
    """

    def __init__(
            self,
            maxfun=15000,
            xtol=-1,
            ftol=-1,
            pgtol=1e-05,
            epsilon=1e-08,
            screen_output=False):
        """
        Constructs a Truncated Newton Method algorithm (SciPy)

        NOTE: gradient is numerically approximated

        USAGE: algorithm.scipy_tnc(maxfun = 1, xtol = -1, ftol = -1, pgtol = 1e-05, epsilon = 1e-08, screen_output = False)

        * maxfun: Maximum number of function evaluation.
        * xtol: Precision goal for the value of x in the stopping criterion
                (after applying x scaling factors). If xtol < 0.0, xtol is set to
                sqrt(machine_precision). Defaults to -1.
        * ftol: Precision goal for the value of f in the stoping criterion.
                If ftol < 0.0, ftol is set to 0.0 defaults to -1.
        * pgtol: Precision goal for the value of the projected gradient in the
                 stopping criterion (after applying x scaling factors). If pgtol
                 < 0.0, pgtol is set to 1e-2 * sqrt(accuracy).
                 Setting it to 0.0 is not recommended. Defaults to -1.
        * epsilon: The stepsize in a finite difference approximation for the objfun
        * screen_output: Set to True to print iterations
        """
        _scipy_base.__init__(self, 'fmin_tnc', False)
        self.maxfun = maxfun
        self.xtol = xtol
        self.ftol = ftol
        self.pgtol = pgtol
        self.epsilon = epsilon
        self.screen_output = screen_output

    def evolve(self, pop):
        from numpy import concatenate, array
        prob = pop.problem
        self._problem_checks(prob)
        if len(pop) == 0:
            return pop
        _, x0, x0_comb = self._starting_params(pop)
        # Extract the a list of tuples representing the bounds.
        prob_bounds = [(prob.lb[i], prob.ub[i])
                       for i in range(0, prob.dimension - prob.i_dimension)]
        if self.screen_output:
            msg = 15
        else:
            msg = 0
        retval = self.solver(
            lambda x: array(
                prob.objfun(
                    concatenate(
                        (x,
                         x0_comb))),
                dtype=float),
            x0,
            bounds=prob_bounds,
            approx_grad=True,
            messages=msg,
            maxfun=self.maxfun,
            xtol=self.xtol,
            ftol=self.ftol,
            pgtol=self.pgtol,
            epsilon=self.epsilon)
        new_chromosome = list(retval[0]) + list(x0_comb)
        pop.set_x(
            pop.get_best_idx(),
            self._check_new_chromosome(
                new_chromosome,
                prob))
        return pop

    def get_name(self):
        return "Truncated Newton Method (SciPy)"

    def human_readable_extra(self):
        return "maxfun = " + str(self.maxfun) + ", xtol = " + str(self.xtol) + ", ftol = " + \
            str(self.ftol) + ", pgtol = " + str(self.pgtol) + ", epsilon = " + str(self.epsilon)


class scipy_cobyla(_scipy_base):

    """
    Wrapper around SciPy's cobyla optimiser.
    """

    def __init__(self, max_fun=1, rho_end=1E-5, screen_output=False):
        """
        Constructs a Constrained Optimization BY Linear Approximation (COBYLA) algorithm (SciPy)

        NOTE: equality constraints are transformed into two inequality constraints automatically

        USAGE: algorithm.scipy_cobyla(max_fun = 1,rho_end = 1E-5,screen_output = False)

        * maxfun: Maximum number of function evaluations.
        * rhoend: Final accuracy in the optimization (not precisely guaranteed). This is a lower bound on the size of the trust region.
        * screen_output: Set to True to print iterations
        """
        _scipy_base.__init__(self, 'fmin_cobyla', True)
        self.max_fun = max_fun
        self.rho_end = rho_end
        self.screen_output = screen_output

    def evolve(self, pop):
        from numpy import concatenate, array
        prob = pop.problem
        self._problem_checks(prob)
        if len(pop) == 0:
            return pop
        nec, x0, x0_comb = self._starting_params(pop)
        # We need to build a vector of functions for the constraints. COBYLA
        # uses inequality constraints with >= 0.
        ub = prob.ub
        lb = prob.lb
        f_cons = []
        for i in range(0, nec):
            # Each eq. constraint is converted into two ineq. constraints with
            # different signs,
            f_cons.append(
                lambda x, i=i: prob.compute_constraints(
                    concatenate(
                        (x, x0_comb)))[i])
            f_cons.append(
                lambda x,
                i=i: -
                prob.compute_constraints(
                    concatenate(
                        (x,
                         x0_comb)))[i])
        for i in range(nec, prob.c_dimension):
            # Ineq. constraints.
            f_cons.append(
                lambda x,
                i=i: -
                prob.compute_constraints(
                    concatenate(
                        (x,
                         x0_comb)))[i])
        for i in range(0, prob.dimension - prob.i_dimension):
            # Box bounds implemented as inequality constraints.
            f_cons.append(lambda x, i=i: x[i] - lb[i])
            f_cons.append(lambda x, i=i: ub[i] - x[i])
        if self.screen_output:
            iprn = 1
        else:
            iprn = 0
        retval = self.solver(
            lambda x: array(
                prob.objfun(
                    concatenate(
                        (x,
                         x0_comb))),
                dtype=float),
            x0,
            cons=f_cons,
            iprint=iprn,
            maxfun=self.maxfun,
            rhoend=self.rhoend)
        new_chromosome = list(retval) + list(x0_comb)
        pop.set_x(
            pop.get_best_idx(),
            self._check_new_chromosome(
                new_chromosome,
                prob))
        return pop

    def get_name(self):
        return "Constrained Optimization BY Linear Approximation (SciPy)"

    def human_readable_extra(self):
        return "maxfun = " + \
            str(self.maxfun) + ", rhoend = " + str(self.rhoend)


# This algorithm suck at the moment and thus I will not include it


# class scipy_anneal(_scipy_base):
#	"""
#	Wrapper around SciPy's anneal optimiser.
#	"""
# def __init__(self, schedule = 'fast', screen_output = False, T0 = None, Tf = 9.9999999999999998e-13, maxfun = 1000000, maxaccept = None, maxiter = 100000, boltzmann #= 1.0, learn_rate = 0.5, feps = 9.9999999999999995e-07, dwell = 5):
#		"""
#		Constructs a Simulate Annealing algorithm
#
#		USAGE: algorithm.scipy_anneal(maxiter = 100000, T0 = None, Tf = 9.9999999999999998e-13, schedule = 'fast',
#		                              maxfun = 1000000, maxaccept = None, boltzmann = 1.0, learn_rate = 0.5,
#		                              feps = 9.9999999999999995e-07, dwell = 50, screen_output = False)
#
#		NOTE: it is not guaranteed that the objective function will be only called with
#		      chrmosomes within the bounds
#
#		* maxiter: Maximum cooling iterations
#		* T0: Starting temperature
#		* Tf: End Temperature
#		* schedule: Annealing schedule one of 'fast', 'cauchy', 'boltzmann'
#		* maxfun: Maximum function evaluations
#		* maxaccept: Maximum change to accept
#		* boltzmann: Boltzmann constant in acceptance test
#                             (increase for less stringent test at each temperature).
#		* learn_rate: Scale constant for adjusting guesses
#		* feps: Stopping relative error tolerance for the function value in
#                        last four coolings.
#		* dwell: The number of times to search the space at each temperature.
#		* screen_output: Set to True to print iterations
#		"""
#
#		_scipy_base.__init__(self,'anneal',False)
#		self.schedule = schedule
#		self.T0 = T0
#		self.Tf = Tf
#		self.maxfun = maxfun
#		self.maxaccept = maxaccept
#		self.maxiter = maxiter
#		self.boltzmann = boltzmann
#		self.learn_rate = learn_rate
#		self.feps = feps
#		self.dwell = dwell
#	def evolve(self,pop):
#		from numpy import concatenate, array
#		prob = pop.problem
#		self._problem_checks(prob)
# If population is empty, just return input population.
#		if len(pop) == 0:
#			return pop
# Get starting params.
#		n_ec, x0, x0_comb = self._starting_params(pop)
# Run the optimisation.
#		retval = self.solver(lambda x: prob.objfun(concatenate((x, x0_comb)))[0],x0,lower = array(prob.lb,dtype=float),upper = array(prob.ub,dtype=float)
#			          , full_output = int(self.screen_output), schedule = self.schedule, T0 = self.T0, Tf = self.Tf, maxeval = self.maxfun
#                                 , maxaccept = self.maxaccept, maxiter = self.maxiter, boltzmann = self.boltzmann, learn_rate = self.learn_rate
#                                  , feps = self.feps, dwell = self.dwell)
# Set the individual's chromosome in the population and return. Conserve the integer part from the
# original individual.
#		new_chromosome = list(retval[0]) + list(x0_comb)
#		pop.set_x(0,self._check_new_chromosome(new_chromosome,prob))
#		return pop
#	def get_name(self):
#		return "Simulated Annealing (SciPy)"
#	def human_readable_extra(self):
#		return ("maxiter = " + str(self.maxiter) + ", T0 = " + str(self.T0) + ", Tf = " + str(self.Tf) + ", schedule = " + self.schedule
#		                    + ", maxfun = " + str(self.maxfun) + ", maxaccept = " + str(self.maxaccept) + ", boltzmann = " + str(self.boltzmann)
#		                    + ", learn_rate = " + str(self.learn_rate) + ", feps = " + str(self.feps) + ", dwell = " + str(self.dwell)
