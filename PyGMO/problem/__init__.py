# -*- coding: iso-8859-1 -*-
from PyGMO import __extensions__
from PyGMO.problem._base import base
from PyGMO.problem._base_stochastic import base_stochastic
from PyGMO.problem._problem import *
from PyGMO.problem._problem_meta import *
try:
    from PyGMO.problem._problem_space import *
    from PyGMO.problem._gtop import *
    __extensions__["gtop"] = True
except ImportError:
    pass
from PyGMO.problem._problem import _base
from PyGMO.problem._problem import _base_stochastic
from PyGMO.problem._example import py_example
from PyGMO.problem._example import py_example_max
from PyGMO.problem._example_stochastic import py_example_stochastic
from PyGMO.problem._pl2pl import py_pl2pl
from PyGMO.problem._mo import *
from PyGMO.problem._tsp import *


# If GSL support is active import mit_sphere
try:
    from PyGMO.problem._mit_spheres import visualize as _visualize
    mit_spheres.visualize = _visualize

    def _mit_spheres_ctor(
        self,
        sample_size=10,
        n_hidden=10,
        ode_prec=1E-3,
        seed=0,
        symmetric=False,
        simulation_time=50.0,
        sides=[
            0.6,
            0.7,
            0.8]):
        """
        Construct a Neurocontroller Evolution problem that seeks to drive three point masses to form a triangle
        This problem was used to design a contorller for the MIT SPHERES test bed on boear the ISS

        USAGE: problem.mit_spheres(sample_size = 10, n_hidden = 10, ode_prec = 1E-3, seed = 0, symmetric = False, simulation_time = 50.0):

        * sample_size: number of initial conditions the neurocontroller is tested from
        * n_hidden: number of hidden  for the feed-forward neural network
        * ode_prec: relative numerical precision of neurons the ODE integrator
        * seed: integer used as starting random seed to build the pseudorandom sequences used to generate the sample
        * symmetric: when True activates a Neural Network having symmetric weights (i.e. purely homogeneuos agents)
        * simulation_time: when True activates a Neural Network having symmetric weights (i.e. purely homogeneuos agents)
        * sides: sides of the triangle

"""

        # We construct the arg list for the original constructor exposed by
        # boost_python
        arg_list = []
        arg_list.append(sample_size)
        arg_list.append(n_hidden)
        arg_list.append(ode_prec)
        arg_list.append(seed)
        arg_list.append(symmetric)
        arg_list.append(simulation_time)
        arg_list.append(sides)
        self._orig_init(*arg_list)
    mit_spheres._orig_init = mit_spheres.__init__
    mit_spheres.__init__ = _mit_spheres_ctor

    from PyGMO import __version__
    __version__ = __version__ + "GTOP " + "GSL "

except:
    pass

# Creating the list of problems


def _get_problem_list():
    from PyGMO import problem
    return [
        problem.__dict__[n] for n in [
            n for n in dir(problem) if not n.startswith('_') and not n == 'base' and not n == "base_stochastic" and (
                issubclass(
                    problem.__dict__[n],
                    problem._base) or issubclass(
                    problem.__dict__[n],
                    problem._base_stochastic))]]

# Redefining the constructors of all problems to obtain good documentation
# and allowing kwargs


def _rastrigin_ctor(self, dim=10):
    """
    Constructs a Rastrigin problem (Box-Constrained Continuous Single-Objective)

    USAGE: problem.rastrigin(dim=10)

    * dim: problem dimension
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(dim)
    self._orig_init(*arg_list)
rastrigin._orig_init = rastrigin.__init__
rastrigin.__init__ = _rastrigin_ctor


def _rosenbrock_ctor(self, dim=10):
    """
    Constructs a Rosenbrock problem (Box-Constrained Continuous Single-Objective)

    USAGE: problem.rosenbrock(dim=10)

    * dim: problem dimension
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(dim)
    self._orig_init(*arg_list)
rosenbrock._orig_init = rosenbrock.__init__
rosenbrock.__init__ = _rosenbrock_ctor


def _ackley_ctor(self, dim=10):
    """
    Constructs a Ackley problem (Box-Constrained Continuous Single-Objective)

    USAGE: problem.ackley(dim=10)

    * dim: problem dimension
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(dim)
    self._orig_init(*arg_list)
ackley._orig_init = ackley.__init__
ackley.__init__ = _ackley_ctor


def _schwefel_ctor(self, dim=10):
    """
    Constructs a Schwefel problem (Box-Constrained Continuous Single-Objective)

    USAGE: problem.schwefel(dim=10)

    * dim: problem dimension
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(dim)
    self._orig_init(*arg_list)
schwefel._orig_init = schwefel.__init__
schwefel.__init__ = _schwefel_ctor


def _dejong_ctor(self, dim=10):
    """
    Constructs a De Jong problem (Box-Constrained Continuous Single-Objective)

    USAGE: problem.dejong(dim=10)

    * dim: problem dimension
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(dim)
    self._orig_init(*arg_list)
dejong._orig_init = dejong.__init__
dejong.__init__ = _dejong_ctor


def _griewank_ctor(self, dim=10):
    """
    Constructs a Griewank problem (Box-Constrained Continuous Single-Objective)

    USAGE: problem.griewank(dim=10)

    * dim: problem dimension
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(dim)
    self._orig_init(*arg_list)
griewank._orig_init = griewank.__init__
griewank.__init__ = _dejong_ctor


def _lavor_maculan_ctor(self, n_atoms=4):
    """
    Constructs a Lavor-Maculan problem (Box-Constrained Continuous Single-Objective) with hydrocarbon chain of N atoms (N-3 dimensions).

    USAGE: problem.lavor_maculan(n_atoms=10)

    * n_atoms: number of atoms
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(n_atoms)
    self._orig_init(*arg_list)
lavor_maculan._orig_init = lavor_maculan.__init__
lavor_maculan.__init__ = _lavor_maculan_ctor


def _lennard_jones_ctor(self, n_atoms=4):
    """
    Constructs a Lennard-Jones problem (Box-Constrained Continuous Single-Objective)

    USAGE: problem.lennard_jones(n_atoms=4)

    * n_atoms: number of atoms
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(n_atoms)
    self._orig_init(*arg_list)
lennard_jones._orig_init = lennard_jones.__init__
lennard_jones.__init__ = _lennard_jones_ctor


def _branin_ctor(self):
    """
    Constructs a Branin problem (Box-Constrained Continuous Single-Objective)

    USAGE: problem.branin()

    """
    self._orig_init()

branin._orig_init = branin.__init__
branin.__init__ = _branin_ctor


def _himmelblau_ctor(self):
    """
    Constructs a Himmelblau problem (Box-Constrained Continuous Single-Objective)

    USAGE: problem.himmelblau()

    """
    self._orig_init()

himmelblau._orig_init = himmelblau.__init__
himmelblau.__init__ = _himmelblau_ctor


def _bukin_ctor(self):
    """
    Constructs a Bukin's f6 problem (Box-Constrained Continuous Single-Objective)

    USAGE: problem.bukin()

    """
    self._orig_init()

bukin._orig_init = bukin.__init__
bukin.__init__ = _bukin_ctor


def _michalewicz_ctor(self, dim=10):
    """
    Constructs a Michalewicz problem (Box-Constrained Continuous Single-Objective)

    USAGE: problem.michalewicz(dim=5)

    NOTE: Minimum is -4.687 for dim=5 and -9.66 for dim = 10

    * dim: problem dimension
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(dim)
    self._orig_init(*arg_list)
michalewicz._orig_init = michalewicz.__init__
michalewicz.__init__ = _michalewicz_ctor


def _kur_ctor(self, dim=10):
    """
    Constructs a Kursawe's study problem (Box-Constrained Continuous Multi-Objective)

    NOTE: K Deb, A Pratap, S Agarwal: A fast and elitist multiobjective genetic algorithm: NSGA-II, IEEE Transactions on, 2002

    USAGE: problem.kur(dim = 10)

    * dim: problem dimension
    """
    arg_list = []
    arg_list.append(dim)
    self._orig_init(*arg_list)

kur._orig_init = kur.__init__
kur.__init__ = _kur_ctor


def _fon_ctor(self):
    """
    Constructs a Fonseca and Fleming's study problem (Box-Constrained Continuous Multi-Objective)

    NOTE: K Deb, A Pratap, S Agarwal: A fast and elitist multiobjective genetic algorithm: NSGA-II, IEEE Transactions on, 2002

    USAGE: problem.fon()
    """
    arg_list = []
    self._orig_init(*arg_list)

fon._orig_init = fon.__init__
fon.__init__ = _fon_ctor


def _pol_ctor(self):
    """
    Constructs a Poloni's study study problem (Box-Constrained Continuous Multi-Objective)

    NOTE: K Deb, A Pratap, S Agarwal: A fast and elitist multiobjective genetic algorithm: NSGA-II, IEEE Transactions on, 2002

    USAGE: problem.pol()
    """
    arg_list = []
    self._orig_init(*arg_list)

pol._orig_init = pol.__init__
pol.__init__ = _pol_ctor


def _sch_ctor(self):
    """
    Constructs a Schaffer's study problem (Box-Constrained Continuous Multi-Objective)

    NOTE: K Deb, A Pratap, S Agarwal: A fast and elitist multiobjective genetic algorithm: NSGA-II, IEEE Transactions on, 2002

    USAGE: problem.sch()
    """
    arg_list = []
    self._orig_init(*arg_list)

sch._orig_init = sch.__init__
sch.__init__ = _sch_ctor


def _pressure_vessel_ctor(self):
    """
    Constructs a pressure vessel design problem (Constrained Continuous Single-Objective)

    USAGE: problem.pressure_vessel()
    """
    arg_list = []
    self._orig_init(*arg_list)

pressure_vessel._orig_init = pressure_vessel.__init__
pressure_vessel.__init__ = _pressure_vessel_ctor


def _tens_comp_string_ctor(self):
    """
    Constructs a tension compression string design problem (Constrained Continuous Single-Objective)

    USAGE: problem.tens_comp_string()
    """
    arg_list = []
    self._orig_init(*arg_list)

tens_comp_string._orig_init = tens_comp_string.__init__
tens_comp_string.__init__ = _tens_comp_string_ctor


def _welded_beam_ctor(self):
    """
    Constructs a welded beam design problem (Constrained Continuous Single-Objective)

    USAGE: problem.welded_beam()
    """
    arg_list = []
    self._orig_init(*arg_list)

welded_beam._orig_init = welded_beam.__init__
welded_beam.__init__ = _welded_beam_ctor


def _cec2006_ctor(self, prob_id=1):
    """
    Constructs one of the 24 CEC2006 Competition Problems (Constrained Continuous Single-Objective)

    USAGE: problem.cec2006(prob_id=1)

    * prob_id: Problem number, one of [1,2,...24]
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(prob_id)
    self._orig_init(*arg_list)

cec2006._orig_init = cec2006.__init__
cec2006.__init__ = _cec2006_ctor


def _cec2009_ctor(self, prob_id=1, dim=30, is_constrained=False):
    """
    Constructs one of the 20 CEC2009 Competition Problems (Constrained / Unconstrained Multi-Objective)

    USAGE: problem.cec2009(prob_id=1, dim=30, is_constrained=False)

    * prob_id: Problem number, one of [1,2,...10]
    * dim: Problem's dimension (default is 30, corresponding to the competition set-up)
    * is_constrained: if True constructs the CF problems, otherwise the UF (constrained/unconstrained)
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(prob_id)
    arg_list.append(dim)
    arg_list.append(is_constrained)
    self._orig_init(*arg_list)

cec2009._orig_init = cec2009.__init__
cec2009.__init__ = _cec2009_ctor


def _cec2013_ctor(self, prob_id=1, dim=10, path="input_data/"):
    """
    Constructs one of the 28 CEC2013 Competition Problems (Box-Constrained Continuous Single-Objective)

    NOTE: this problem requires two files to be put in the path indicated: "M_Dxx.txt" and "shift_data.txt".
    These files can be downloaded from the CEC2013 competition site: http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/CEC2013/cec13-c-code.zip

    USAGE: problem.cec2013(dim = 10, prob_id=1, path="input_data/")

    * prob_id: Problem number, one of [1,2,...10]
    * dim: Problem's dimension (default is 10)
    * path: Whether the problem is constrained or unconstrained

    """
    arg_list = []
    arg_list.append(prob_id)
    arg_list.append(dim)
    arg_list.append(path)
    self._orig_init(*arg_list)

cec2013._orig_init = cec2013.__init__
cec2013.__init__ = _cec2013_ctor


def _luksan_vlcek_1_ctor(self, dim=3):
    """
    Constructs the first Luksan Vlcek problem (Constrained Continuous Single-Objective)

    NOTE: L. Luksan and J. Vlcek, "Sparse and Parially Separable Test Problems for Unconstrained and Equality Constrained Optimization"

    USAGE: problem.luksan_vlcek_1(dim=3)

    * dim: problem dimension
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(dim)
    self._orig_init(*arg_list)
luksan_vlcek_1._orig_init = luksan_vlcek_1.__init__
luksan_vlcek_1.__init__ = _luksan_vlcek_1_ctor


def _luksan_vlcek_2_ctor(self, dim=16):
    """
    Constructs the second Luksan Vlcek problem (Constrained Continuous Single-Objective)

    NOTE: L. Luksan and J. Vlcek, "Sparse and Parially Separable Test Problems for Unconstrained and Equality Constrained Optimization"

    USAGE: problem.luksan_vlcek_2(dim=16)

    * dim: problem dimension
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(dim)
    self._orig_init(*arg_list)
luksan_vlcek_2._orig_init = luksan_vlcek_2.__init__
luksan_vlcek_2.__init__ = _luksan_vlcek_2_ctor


def _luksan_vlcek_3_ctor(self, dim=8):
    """
    Constructs the third Luksan Vlcek problem (Constrained Continuous Single-Objective)

    NOTE: L. Luksan and J. Vlcek, "Sparse and Parially Separable Test Problems for Unconstrained and Equality Constrained Optimization"

    USAGE: problem.luksan_vlcek_3(dim=8)

    * dim: problem dimension
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(dim)
    self._orig_init(*arg_list)
luksan_vlcek_3._orig_init = luksan_vlcek_3.__init__
luksan_vlcek_3.__init__ = _luksan_vlcek_3_ctor


def _snopt_toyprob_ctor(self):
    """
    Constructs SNOPT toy-problem (Box-Constrained Continuous Multi-Objective)

    USAGE: problem.snopt_toyprob()
    """
    arg_list = []
    self._orig_init(*arg_list)

snopt_toyprob._orig_init = snopt_toyprob.__init__
snopt_toyprob.__init__ = _snopt_toyprob_ctor


def _string_match_ctor(self, string="Can we use it for space?"):
    """
    Constructs a string-match problem (Box-Constrained Integer Single-Objective)

    NOTE: This is the problem of matching a string. Transcribed as an optimization problem

    USAGE: problem.string_match(string = "mah")

    * string: string to match
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(string)
    self._orig_init(*arg_list)
string_match._orig_init = string_match.__init__
string_match.__init__ = _string_match_ctor


def _golomb_ruler_ctor(self, order=5, length=10):
    """
    Constructs a Golomb Ruler problem (Constrained Integer Single-Objective)

    NOTE: see http://en.wikipedia.org/wiki/Golomb_ruler

    USAGE: problem.golomb_ruler(order = 5, length=10)

    * order: order of the Golomb ruler
    * length: length of the Golomb ruler
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(order)
    arg_list.append(length)
    self._orig_init(*arg_list)
golomb_ruler._orig_init = golomb_ruler.__init__
golomb_ruler.__init__ = _golomb_ruler_ctor


def _inventory_ctor(self, weeks=4, sample_size=10, seed=0):
    """
    Constructs an Inventory Problem (Stochastic Objective Function)

    NOTE: see www2.isye.gatech.edu/people/faculty/Alex_Shapiro/SPbook.pdf

    USAGE: problem.inventory(weeks = 4, sample_size = 10, seed = 0):

    * week: dimension of the problem corresponding to the numer of weeks
                     to plan the inventory for.
    * sample_size: dimension of the sample used to approximate the expected value
    * seed: integer used as starting random seed to build the
                     pseudorandom sequences used to generate the sample
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(weeks)
    arg_list.append(sample_size)
    arg_list.append(seed)
    self._orig_init(*arg_list)
inventory._orig_init = inventory.__init__
inventory.__init__ = _inventory_ctor


def _normalized_ctor(self, problem=None):
    """
    Normalizes a problem (e.g. maps all variables to [-1,1])

    NOTE: this meta-problem constructs a new problem having normalized bounds/variables

    USAGE: problem.normalized(problem=PyGMO.ackley(1))

    * problem: PyGMO problem one wants to normalize

    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    if problem is None:
        problem = ackley(1)
    arg_list.append(problem)
    self._orig_init(*arg_list)
normalized._orig_init = normalized.__init__
normalized.__init__ = _normalized_ctor


def _decompose_ctor(
        self,
        problem=None,
        method='tchebycheff',
        weights=[],
        z=[]):
    """
    Implements a meta-problem class resulting in a decomposed version
    of the multi-objective input problem, i.e. a single-objective problem
    having as fitness function some kind of combination of the original fitness functions.

    NOTE: this meta-problem constructs a new single-objective problem

    USAGE: problem.decompose(problem=PyGMO.zdt(1, 2), method = 'tchebycheff', weights=a random vector (summing to one), z= a zero vector)

    * problem: PyGMO problem one wants to decompose
    * method: the decomposition method to use ('weighted', 'tchebycheff' or 'bi')
    * weights: the weight vector to build the new fitness function
    * z: the reference point (used in TCHEBYCHEFF and BI methods)

    """

    # We construct the arg list for the original constructor exposed by
    # boost_python

    DECOMPOSITION_TYPE = {
        'weighted': _problem_meta._decomposition_method.WEIGHTED,
        'tchebycheff': _problem_meta._decomposition_method.TCHEBYCHEFF,
        'bi': _problem_meta._decomposition_method.BI,
    }

    arg_list = []
    if problem is None:
        problem = zdt(1, 2)
    arg_list.append(problem)
    arg_list.append(DECOMPOSITION_TYPE[method.lower()])
    arg_list.append(weights)
    arg_list.append(z)
    self._orig_init(*arg_list)
decompose._orig_init = decompose.__init__
decompose.__init__ = _decompose_ctor


def _shifted_ctor(self, problem=None, shift=None):
    """
    Shifts a problem.

    NOTE: this meta-problem constructs a new problem where the objective function will be f(x+b),
          where b is the shift (bounds are also chaged accordingly)

    USAGE: problem.shifted(problem=PyGMO.ackley(1), shift = a random vector)

    * problem: PyGMO problem one wants to shift
    * shift: a value or a list containing the shifts. By default, a radnom shift is created within the problem bounds

    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    if problem is None:
        problem = ackley(1)
    arg_list.append(problem)
    if shift is not None:
        arg_list.append(shift)
    self._orig_init(*arg_list)
shifted._orig_init = shifted.__init__
shifted.__init__ = _shifted_ctor


def _rotated_ctor(self, problem=None, rotation=None):
    """
    Rotates a problem. (also reflections are possible)
    The new objective function will be f(Rx_{normal}), where R is an orthogonal matrix and x_{normal}
    is the decision vector normailized to [-1,1]

    NOTE: To ensure all of the original space is included in the new box-constrained search space, bounds
    of the normalized variables are expanded to [-sqrt(2),sqrt(2)]. It is still guaranteed theat the original
    objective function will not be called outside of the original bounds by projecting points outside the original
    space onto the boundary

    USAGE: problem.rotated(problem=PyGMO.ackley(1), rotation = a random orthogonal matrix)

    * problem: PyGMO problem one wants to rotate
    * rotation: a list of lists (matrix). If not specified, a random orthogonal matrix is used.

    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    if problem is None:
        problem = ackley(1)
    arg_list.append(problem)
    if rotation is not None:
        arg_list.append(rotation)
    self._orig_init(*arg_list)
rotated._orig_init = rotated.__init__
rotated.__init__ = _rotated_ctor


_problem_meta.noisy.noise_distribution = _problem_meta._noise_distribution


def _noisy_ctor(
        self,
        problem=None,
        trials=1,
        param_first=0.0,
        param_second=1.0,
        noise_type=noisy.noise_distribution.NORMAL,
        seed=0):
    """
    Inject noise to a problem.
    The new objective function will become stochastic, influence by a normally distributed noise.

    USAGE: problem.noisy(problem=PyGMO.ackley(1), trials = 1, param_first=0.0, param_second=1.0, noise_type = problem.noisy.noise_distribution.NORMAL, seed=0)

    * problem: PyGMO problem on which one wants to add noises
    * trials: number of trials to average around
    * param_first: Mean of the Gaussian noise / Lower bound of the uniform noise
    * param_second: Standard deviation of the Gaussian noise / Upper bound of the uniform noise
    * noise_type: Whether to inject a normally distributed noise or uniformly distributed noise
    * seed: Seed for the underlying RNG
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    if problem is None:
        problem = ackley(1)
    arg_list.append(problem)
    arg_list.append(trials)
    arg_list.append(param_first)
    arg_list.append(param_second)
    arg_list.append(noise_type)
    arg_list.append(seed)
    self._orig_init(*arg_list)
noisy._orig_init = noisy.__init__
noisy.__init__ = _noisy_ctor


def _robust_ctor(self, problem=None, trials=1, rho=0.1, seed=0):
    """
    Inject noise to a problem in the decision space.
    The solution to the resulting problem is robust the the noise in the rho area.

    USAGE: problem.robust(problem=PyGMO.ackley(10), trials=1, rho=0.1, seed=0)

    * problem: PyGMO problem to be transformed to its robust version
    * trials: number of trials to average around
    * rho: Parameter controlling the magnitude of noise
    * seed: Seed for the underlying RNG
    """
    arg_list = []
    if problem is None:
        problem = ackley(10)
    arg_list.append(problem)
    arg_list.append(trials)
    arg_list.append(rho)
    arg_list.append(seed)
    self._orig_init(*arg_list)
robust._orig_init = robust.__init__
robust.__init__ = _robust_ctor

# Renaming and placing the enums
_problem_meta.death_penalty.method = _problem_meta._death_method_type


def _death_penalty_ctor(self, problem=None, method=None, penalty_factors=None):
    """
    Implements a meta-problem class that wraps some other constrained problems, resulting in death penalty constraints handling.
    Three implementations of the death penalty are available. The first one is the most common simple death penalty. The second one is the death
    penalty defined by Angel Kuri Morales et al. (Kuri Morales, A. and Quezada, C.C. A Universal eclectic genetic algorithm for constrained optimization, Proceedings 6th European Congress on Intelligent Techniques & Soft Computing, EUFIT'98, 518-522, 1998.)
    Simple death penalty penalizes the fitness function with a high value, Kuri method penalizes the
    fitness function according to the rate of satisfied constraints. The third one is a weighted static penalization.
    It penalizes the objective with the sum of the constraints violation, each one penalized with a given factor.

    USAGE: problem.death_penalty(problem=PyGMO.cec2006(4), method=death_penalty.method.SIMPLE)

    * problem: PyGMO constrained problem one wants to treat with a death penalty approach
    * method: Simple death method set with SIMPLE, Kuri method set with KURI, weighted static penalization with WEIGHTED
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    if problem is None:
        problem = cec2006(4)
    if method is None:
        method = death_penalty.method.SIMPLE
    arg_list.append(problem)
    arg_list.append(method)
    if penalty_factors is not None:
        arg_list.append(penalty_factors)
    self._orig_init(*arg_list)
death_penalty._orig_init = death_penalty.__init__
death_penalty.__init__ = _death_penalty_ctor


# Renaming and placing the enums
_problem_meta.con2mo.method = _problem_meta._con2mo_method_type


def _con2mo_ctor(self, problem=None, method='obj_cstrsvio'):
    """
    Transforms a constrained problem into a multi-objective problem

    Three implementations of the constrained to multi-objective are available.
    1) 'obj_cstrs': The multi-objective problem is created with two objectives. The first
    objective is the same as that of the input problem, the second is the number of constraint violated
    2) 'obj_cstrsvio': The multi-objective problem is created with two objectives. The first
    objective is the same as that of the input problem, the second is the norm of the total constraint violation
    3) 'obj_eqvio_ineqvio': 	2) 'obj_cstrsvio': The multi-objective problem is created with three objectives. The first
    objective is the same as that of the input problem, the second is the norm of the total equality constraint violation,
    the third is the norm of the total inequality constraint violation.

    USAGE: problem.con2mo(problem=PyGMO.cec2006(4), method='obj_cstrsvio')

    * problem: original PyGMO constrained problem
    * method: one of 'obj_cstrsvio', 'obj_eqvio_ineqvio', 'obj_cstrs'
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python

    METHOD_TYPE = {
        'obj_cstrs': _problem_meta._con2mo_method_type.OBJ_CSTRS,
        'obj_cstrsvio': _problem_meta._con2mo_method_type.OBJ_CSTRSVIO,
        'obj_eqvio_ineqvio': _problem_meta._con2mo_method_type.OBJ_EQVIO_INEQVIO,
    }

    arg_list = []
    if problem is None:
        problem = cec2006(4)
    method = METHOD_TYPE[method.lower()]
    arg_list.append(problem)
    arg_list.append(method)
    self._orig_init(*arg_list)

con2mo._orig_init = con2mo.__init__
con2mo.__init__ = _con2mo_ctor


def _con2uncon_ctor(self, problem=cec2006(4), method='optimality'):
    """
    Implements a meta-problem class that wraps constrained problems,
    resulting in an unconstrained problem. Two methods
    are available for definig the objective function of the meta-problem: 'optimality' and 'feasibility'.
    The 'optimality' uses as objective function the original objective function, it basically removes the constraints from the original problem. The
    'feasibility' uses as objective function the sum of the violation of the constraints, the meta-problem hence optimize just the level of infeasibility.

    Implements a meta-problem class that wraps some other constrained problems,
    resulting in multi-objective problem.

    USAGE: problem.con2uncon(problem=PyGMO.cec2006(4), method='optimality')

    * problem: original PyGMO constrained problem
    * method: one of 'optimality', 'feasibility'.
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    METHOD_TYPE = {
        'optimality': _problem_meta._con2uncon_method_type.OPTIMALITY,
        'feasibility': _problem_meta._con2uncon_method_type.FEASIBILITY,
    }

    arg_list = []
    arg_list.append(problem)
    arg_list.append(METHOD_TYPE[method.lower()])
    self._orig_init(*arg_list)

con2uncon._orig_init = con2uncon.__init__
con2uncon.__init__ = _con2uncon_ctor
