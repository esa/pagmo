# -*- coding: utf-8 -*-
from PyGMO.util._util import *
from PyGMO.util._analysis import *
from PyGMO.util._tsp import read_tsplib

__all__ = ['hypervolume', 'hv_algorithm', 'tsp']


hv_algorithm.__doc__ = """Module containing available algorithms for the hypervolume computation

        USAGE:
                hv_algorithm.hv2d()
                hv_algorithm.hv3d()
                hv_algorithm.hv4d()
                hv_algorithm.wfg()
                hv_algorithm.bf_approx()
                hv_algorithm.bf_fpras()
                hv_algorithm.hoy()
                hv_algorithm.fpl()
"""


class _HypervolumeValidation:

    """
    Utility class containing commonly raised errors.
    Kept in one place to simplify the consistency of error messages across methods
    """

    # Raised when the reference point type is not a list or a tuple, e.g. r =
    # "Foo"
    err_rp_type = TypeError(
        "Reference point must be a list/tuple of real numbers, e.g.: r = [1.0, 1.0, 1.0]")

    # Raised when the reference point is a tuple/list but the items are
    # non-castable to float, e.g. r = [1.0, 2.0, 'foo']
    err_rp_items_type = TypeError(
        "Every item in reference point list/tuple must be castable to float, e.g.: r = [1, '2.5', 10e-4]")

    # Raised when the user does not provide a reference point (mandatory in
    # for every method)
    err_rp_none = TypeError(
        "Reference point (keyword argument 'r') is mandatory")

    # Raised when the user provides something weird as a hv_algorithm, e.g.
    # hv.compute(r=refp, hv_algorithm="A string")
    err_hv_type = TypeError(
        "Hypervolume algorithm must be an instance of a correct type, e.g.: algo = hv_algorithm.wfg()")

    # Raised when the hypervolume object is constructed by anything other than
    # a population object, tuple or a list, e.g. hypervolume("foo bar"),
    # hypervolume([[1,2],[2,"foo"]]) etc.
    err_hv_ctor_type = TypeError(
        "Hypervolume object must be constructed from a list/tuple of points or a population object")

    # Raised when the hypervolume object is constructed with an incorrect
    # keyword argument
    err_hv_ctor_args = TypeError(
        "Hypervolume takes either exactly one unnamed argument or one keyword argument 'data_src' in the constructor")

    # types of hypervolume algorithms
    types_hv_algo = (hv_algorithm.hv2d, hv_algorithm.hv3d, hv_algorithm.hv4d, hv_algorithm.wfg, hv_algorithm.bf_approx, hv_algorithm.bf_fpras, hv_algorithm.hoy, hv_algorithm.fpl)

    # allowed types for the refernce point
    types_rp = (list, tuple,)

    @classmethod
    def handle_refpoint(cls, hypvol, r):
        """
        Common way of validation for the reference point being passed as parameter to 'compute', 'exclusive' and 'least_contributor' methods.
        1. Check if user provided the reference point (mandatory)
        2. Make sure that the reference point is of correct type
        3. Make sure that items of the reference point vector are castable to float
        """
        if r:
            if not any(isinstance(r, T) for T in cls.types_rp):
                raise cls.err_rp_type
            try:
                r = [float(ri) for ri in r]
            except ValueError:
                raise cls.err_rp_items_type
        else:
            raise cls.err_rp_none
        return r

    @classmethod
    def validate_hv_algorithm(cls, algorithm):
        """
        Common way of validation for the hv_algorithm object being passed as parameter to 'compute', 'exclusive' and 'least_contributor' methods, as well as the SMS-EMOA algorithm.
        """
        if not any(isinstance(algorithm, T) for T in cls.types_hv_algo):
            raise cls.err_hv_type
        return algorithm


def _hypervolume_ctor(self, data_src=None, verify=True, *args, **kwargs):
    """
    Constructs a hypervolume object used for the computation of hypervolue and exclusive hypervolume.

    Object can be constructed from the population object, or from a fixed list/tuple of points
    Points within a fixed list must all be of equal size, of dimension larger than 1.

    USAGE:
            from PyGMO import *
            from PyGMO.util import *
            hv = hypervolume(pop)  # Constructs the hypervolume object from individual's fitness vectors
            hv = hypervolume([[1,1,2],[2,1,2],[2,2,3]])
            hv = hypervolume(((1,2), (3,0.5), (1.5, 1.5))
            hv = hypervolume(data_src = ((1,2),(2,3)), verify=False)
    """
    if not data_src or len(args) > 0 or len(kwargs) > 0:
        raise _HypervolumeValidation.err_hv_ctor_args

    from PyGMO import population
    allowed_types = (population, list, tuple,)
    if not any(isinstance(data_src, T) for T in allowed_types):
        raise _HypervolumeValidation.err_hv_ctor_type

    args = []
    args.append(data_src)
    args.append(verify)
    try:
        return self._original_init(*args)
    except TypeError:
        raise _HypervolumeValidation.err_hv_ctor_type
hypervolume._original_init = hypervolume.__init__
hypervolume.__init__ = _hypervolume_ctor


def _hypervolume_compute(self, r=None, algorithm=None, *args, **kwargs):
    """
    Compute the hypervolume indicator for a given reference point, using the provided hypervolume algorithm.
    Type 'hv_algorithm?' for a list of available hypervolume algorithms.

    USAGE:
            hv.compute(r=[5.0]*2)
            hv.compute(r=[5.0]*2, algorithm = hv_algorithm.hv2d())
            * r - reference point used for computation
            * algorithm (optional) - hypervolume algorithm used for the computation, uses the best performing algorithm for given dimension by default.
    """
    if len(args) > 0 or len(kwargs) > 0:
        raise TypeError(
            "Incorrect combination of args/kwargs, type 'hypervolume.compute?' for usage")

    r = _HypervolumeValidation.handle_refpoint(self, r)
    args = []
    args.append(r)
    if algorithm:
        algorithm = _HypervolumeValidation.validate_hv_algorithm(algorithm)
        args.append(algorithm)
    return self._original_compute(*args)

hypervolume._original_compute = hypervolume.compute
hypervolume.compute = _hypervolume_compute


def _hypervolume_exclusive(
        self,
        p_idx=None,
        r=None,
        algorithm=None,
        *args,
        **kwargs):
    """
    Compute the exlusive contribution to the total hypervolume by the point at index p_idx, given a reference point and the provided hypervolume algorithm.
    Type 'hv_algorithm?' for a list of available hypervolume algorithms.

    USAGE:
            hv.exclusive(p_idx=0, r=[5.0]*2)
            hv.exclusive(p_idx=0, r=[5.0]*2, algorithm=hv_algorithm.hv2d())
            * p_idx - index of the point
            * r - reference point used for computation
            * algorithm (optional) - hypervolume algorithm used for the computation, uses the best performing algorithm for given dimension by default
    """
    if p_idx is None:
        raise TypeError(
            "p_idx (non-negative integer) argument is required for computation, type 'hypervolume.exclusive?' for usage.")
    if len(args) > 0 or len(kwargs) > 0:
        raise TypeError(
            "Incorrect combination of args/kwargs, type 'hypervolume.exclusive?' for usage.")

    if not isinstance(p_idx, int) or p_idx < 0:
        raise TypeError(
            "individual index (p_idx) must be a non-negative integer")

    r = _HypervolumeValidation.handle_refpoint(self, r)

    args = []
    args.append(p_idx)
    args.append(r)
    if algorithm:
        algorithm = _HypervolumeValidation.validate_hv_algorithm(algorithm)
        args.append(algorithm)
    return self._original_exclusive(*args)

hypervolume._original_exclusive = hypervolume.exclusive
hypervolume.exclusive = _hypervolume_exclusive


def _hypervolume_least_contributor(
        self,
        r=None,
        algorithm=None,
        *args,
        **kwargs):
    """
    Find the least contributing point among the pareto front approximation.
    Type 'hv_algorithm?' for a list of available hypervolume algorithms.

    USAGE:
            hv.least_contributor(r=[5.0]*3)
            hv.least_contributor(r=[5.0]*3, algorithm=hv_algorithm.hv3d())
            * r - reference point used for computation
            * algorithm (optional) - hypervolume algorithm used for the computation, uses the best performing algorithm for given dimension by default
    """

    if len(args) > 0 or len(kwargs) > 0:
        raise TypeError(
            "Incorrect combination of args/kwargs, type 'hypervolume.least_contributor?' for usage")
    r = _HypervolumeValidation.handle_refpoint(self, r)
    args = []
    args.append(r)
    if algorithm:
        algorithm = _HypervolumeValidation.validate_hv_algorithm(algorithm)
        args.append(algorithm)
    return self._original_least_contributor(*args)

hypervolume._original_least_contributor = hypervolume.least_contributor
hypervolume.least_contributor = _hypervolume_least_contributor


def _hypervolume_greatest_contributor(
        self,
        r=None,
        algorithm=None,
        *args,
        **kwargs):
    """
    Find the least contributing point among the pareto front approximation.
    Type 'hv_algorithm?' for a list of available hypervolume algorithms.

    USAGE:
            hv.greatest_contributor(r=[5.0]*3)
            hv.greatest_contributor(r=[5.0]*3, algorithm=hv_algorithm.hv3d())
            * r - reference point used for computation
            * algorithm (optional) - hypervolume algorithm used for the computation, uses the best performing algorithm for given dimension by default
    """

    if len(args) > 0 or len(kwargs) > 0:
        raise TypeError(
            "Incorrect combination of args/kwargs, type 'hypervolume.greatest_contributor?' for usage")
    r = _HypervolumeValidation.handle_refpoint(self, r)
    args = []
    args.append(r)
    if algorithm:
        algorithm = _HypervolumeValidation.validate_hv_algorithm(algorithm)
        args.append(algorithm)
    return self._original_greatest_contributor(*args)

hypervolume._original_greatest_contributor = hypervolume.greatest_contributor
hypervolume.greatest_contributor = _hypervolume_greatest_contributor


def _hypervolume_contributions(self, r=None, algorithm=None, *args, **kwargs):
    """
    Find the contributions to the hypervolume by each point.
    Type 'hv_algorithm?' for a list of available hypervolume algorithms.

    USAGE:
            hv.contributions(r=[5.0]*3)
            hv.contributions(r=[5.0]*3, algorithm=hv_algorithm.hv3d())
            * r - reference point used for computation
            * algorithm (optional) - hypervolume algorithm used for the computation, uses the best performing algorithm for given dimension by default
    """

    if len(args) > 0 or len(kwargs) > 0:
        raise TypeError(
            "Incorrect combination of args/kwargs, type 'hypervolume.contributions?' for usage")
    r = _HypervolumeValidation.handle_refpoint(self, r)
    args = []
    args.append(r)
    if algorithm:
        algorithm = _HypervolumeValidation.validate_hv_algorithm(algorithm)
        args.append(algorithm)
    return self._original_contributions(*args)

hypervolume._original_contributions = hypervolume.contributions
hypervolume.contributions = _hypervolume_contributions


def _hypervolume_get_nadir_point(self, eps=0.0):
    """
    Return Nadir point for given set of points.

    USAGE:
            hv.nadir_point(eps = 10.0)
            * eps (optional) - value added to every objective in order to assert a strong dominance of reference point (1.0 by default).
    """
    try:
        eps = float(eps)
    except ValueError:
        raise TypeError("Epsilon must be castable to float.")

    if eps < 0.0:
        raise ValueError("Epsilon must be a positive value.")

    return self._original_get_nadir_point(eps)

hypervolume._original_get_nadir_point = hypervolume.get_nadir_point
hypervolume.get_nadir_point = _hypervolume_get_nadir_point


def _hypervolume_set_copy_points(self, b):
    """
    Determine whether the hypervolume object should make a copy of points before doing the computation.
    In most cases, only the first computation behaves as expected.
    Used in cases when the hypervolume object is to be used a 'single use' instance only.

    USAGE:
            hv.set_copy_points(True)
    """
    if not isinstance(b, bool):
        raise TypeError("Argument must be of type 'bool'")

    return self._original_set_copy_points(b)

hypervolume._original_set_copy_points = hypervolume.set_copy_points
hypervolume.set_copy_points = _hypervolume_set_copy_points


def _hypervolume_set_verify(self, b):
    """
    Determines whether the hypervolume object should verify whether the set of points and the reference point meet certain conditions.

    USAGE:
            hv.set_verify(True)
    """
    if not isinstance(b, bool):
        raise TypeError("Argument must be of type 'bool'")

    return self._original_set_verify(b)

hypervolume._original_set_verify = hypervolume.set_verify
hypervolume.set_verify = _hypervolume_set_verify


def _hv2d_ctor(self):
    """
    Hypervolume algorithm: hv2d.
    Computational complexity: O(n*logn)
    Applicable to hypervolume computation problems of dimension=2

    USAGE:
            hv = hypervolume(...) # see 'hypervolume?' for usage
            refpoint=[1.0]*2
            hv.compute(r=refpoint, algorithm=hv_algorithm.hv2d())
            hv.exclusive(p_idx=13, refpoint, algorithm=hv_algorithm.hv2d())
            hv.least_contributor(r=refpoint, algorithm=hv_algorithm.hv2d())
    """
    return self._original_init()
hv_algorithm.hv2d._original_init = hv_algorithm.hv2d.__init__
hv_algorithm.hv2d.__init__ = _hv2d_ctor


def _hv3d_ctor(self):
    """
    Hypervolume algorithm: hv3d.
    This class contains the implementation of efficient hypervolume algorithms for 3 dimensions.
    Computational complexity: O(n*logn)
    Applicable to hypervolume computation problems of dimension=3

    REF: "On the Complexity of Computing the Hypervolume Indicator", Nicola Beume, Carlos M. Fonseca, Manuel Lopez-Ibanez,
    Luis Paquete, Jan Vahrenhold.  IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTTATION. VOL. 13, NO. 5, OCTOBER 2009

    REF: "Computing hypervolume contribution in low dimensions: asymptotically optimal algorithm and complexity results", Michael T. M. Emmerich, Carlos M. Fonseca

    USAGE:
            hv = hypervolume(...) # see 'hypervolume?' for usage
            refpoint = [1.0]*3
            hv.compute(r=refpoint, algorithm=hv_algorithm.hv3d())
            hv.exclusive(p_idx=13, r=refpoint, algorithm=hv_algorithm.hv3d())
            hv.least_contributor(r=refpoint, algorithm=hv_algorithm.hv3d())
    """
    return self._original_init()
hv_algorithm.hv3d._original_init = hv_algorithm.hv3d.__init__
hv_algorithm.hv3d.__init__ = _hv3d_ctor


def _hv4d_ctor(self):
    """
    Hypervolume algorithm: HV4d.
    Computational complexity: O(n^2)
    Applicable to hypervolume computation problems of dimension=4

    REF: Andreia P. Guerreiro, Carlos M. Fonseca, Michael T. Emmerich, "A Fast Dimension-Sweep Algorithm for the Hypervolume Indicator in Four Dimensions",
    CCCG 2012, Charlottetown, P.E.I., August 8–10, 2012.

    USAGE:
            hv = hypervolume(...) # see 'hypervolume?' for usage
            refpoint = [1.0]*4
            hv.compute(r=refpoint, algorithm=hv_algorithm.hv4d())
            hv.exclusive(p_idx=13, r=refpoint, algorithm=hv_algorithm.hv4d())
            hv.least_contributor(r=refpoint, algorithm=hv_algorithm.hv4d())
    """
    return self._original_init()
hv_algorithm.hv4d._original_init = hv_algorithm.hv4d.__init__
hv_algorithm.hv4d.__init__ = _hv4d_ctor


def _fpl_ctor(self):
    """
    Hypervolume algorithm: FPL.
    Computational complexity: O(n ^ (d - 2) * log(n))
    Applicable to hypervolume computation problems of arbitrary dimension.

    USAGE:
            hv = hypervolume(...) # see 'hypervolume?' for usage
            refpoint = [1.0]*5
            hv.compute(r=refpoint, algorithm=hv_algorithm.fpl())
            hv.exclusive(p_idx=13, r=refpoint, algorithm=hv_algorithm.fpl())
            hv.least_contributor(r=refpoint, algorithm=hv_algorithm.fpl())
    """
    return self._original_init()
hv_algorithm.fpl._original_init = hv_algorithm.fpl.__init__
hv_algorithm.fpl.__init__ = _fpl_ctor


def _hoy_ctor(self):
    """
    Hypervolume algorithm: HOY.
    Computational complexity: O(n * log(n) + n ^ (d / 2))
    Applicable to hypervolume computation problems of dimension in [2, ..]


    REF: Nicola Beume and Guenter Rudolph, "Faster S-Metric Calculation by Considering Dominated Hypervolume as Klee's Measure Problem.",
    In: B. Kovalerchuk (ed.): Proceedings of the Second IASTED Conference on Computational Intelligence (CI 2006), pp. 231-236.  ACTA Press: Anaheim, 2006.

    USAGE:
            hv = hypervolume(...) # see 'hypervolume?' for usage
            refpoint = [1.0]*5
            hv.compute(r=refpoint, algorithm=hv_algorithm.hoy())
            hv.exclusive(p_idx=13, r=refpoint, algorithm=hv_algorithm.hoy())
            hv.least_contributor(r=refpoint, algorithm=hv_algorithm.hoy())
    """
    return self._original_init()
hv_algorithm.hoy._original_init = hv_algorithm.hoy.__init__
hv_algorithm.hoy.__init__ = _hoy_ctor


def _wfg_ctor(self, stop_dimension=2):
    """
    Hypervolume algorithm: WFG.
    Applicable to hypervolume computation problems of dimension in [2, ..]

    REF: "A Fast Way of Calculating Exact Hypervolumes", Lyndon While, Lucas Bradstreet, Luigi Barone.
    IEEE TRANSACXTIONS ON EVOLUTIONARY COMPUTATION, VOL. 16, NO. 1, FEBRURARY 2012

    USAGE:
            hv = hypervolume(...) # see 'hypervolume?' for usage
            refpoint = [1.0]*7
            hv.compute(r=refpoint, algorithm=hv_algorithm.wfg())
            hv.exclusive(p_idx=13,r=refpoint, algorithm=hv_algorithm.wfg())
            hv.least_contributor(r=refpoint, algorithm=hv_algorithm.wfg())
    """
    args = []
    args.append(stop_dimension)
    return self._original_init(*args)
hv_algorithm.wfg._original_init = hv_algorithm.wfg.__init__
hv_algorithm.wfg.__init__ = _wfg_ctor


def _bf_approx_ctor(
        self,
        use_exact=True,
        trivial_subcase_size=1,
        eps=1e-2,
        delta=1e-6,
        gamma=0.25,
        delta_multiplier=0.775,
        initial_delta_coeff=0.1,
        alpha=0.2):
    """
    Hypervolume algorithm: Bringmann-Friedrich approximation.

    Default values for the parameters of the algorithm were obtained from the shark implementation of the algorithm:
            http://image.diku.dk/shark/doxygen_pages/html/_least_contributor_approximator_8hpp_source.html

    REF: "Approximating the least hypervolume contributor: NP-hard in general, but fast in practice", Karl Bringmann, Tobias Friedrich.

    USAGE:
            * use_exact - should bf_approx use exact methods for computation
            * trivial_subcase_size - when the number of points overlapping the bounding box is smaller or equal to that argument, we compute the exlusive hypervolume exactly
            * eps - accuracy of approximation
            * delta - confidence of approximation
            * gamma - constant used for computation of delta for each of the points during the sampling
            * delta_multiplier - factor with which delta diminishes each round
            * initial_delta_coeff - initial coefficient multiplied by the delta at round 0
            * alpha - coefficicient stating how accurately current lowest contributor should be sampled
            hv = hypervolume(...) # see 'hypervolume?' for usage
            refpoint = [1.0]*7
            hv.least_contributor(r=refpoint, algorithm=hv_algorithm.bf_approx())
    """

    args = []
    args.append(use_exact)
    args.append(trivial_subcase_size)
    args.append(eps)
    args.append(delta)
    args.append(delta_multiplier)
    args.append(alpha)
    args.append(initial_delta_coeff)
    args.append(gamma)
    return self._original_init(*args)
hv_algorithm.bf_approx._original_init = hv_algorithm.bf_approx.__init__
hv_algorithm.bf_approx.__init__ = _bf_approx_ctor


def _bf_fpras_ctor(self, eps=1e-2, delta=1e-2):
    """
    Hypervolume algorithm: Bringmann-Friedrich approximation.

    Default values for the parameters of the algorithm were obtained from the shark implementation of the algorithm:
            http://image.diku.dk/shark/doxygen_pages/html/_least_contributor_approximator_8hpp_source.html

    REF: "Approximating the volume of unions and intersections of high-dimensional geometric objects", Karl Bringmann, Tobias Friedrich.

    USAGE:
            * eps - accuracy of approximation
            * delta - confidence of approximation

            hv = hypervolume(...) # see 'hypervolume?' for usage
            refpoint = [1.0]*7
            hv.least_contributor(r=refpoint, algorithm=hv_algorithm.bf_fpras())
    """

    args = []
    args.append(eps)
    args.append(delta)
    return self._original_init(*args)
hv_algorithm.bf_fpras._original_init = hv_algorithm.bf_fpras.__init__
hv_algorithm.bf_fpras.__init__ = _bf_fpras_ctor


def _race_pop_ctor(self, pop=None, seed=0):
    """
    Constructs a racing object responsible for racing individuals in a population

    USAGE: race_pop(pop, seed=0)

    * pop: The pop containing the individuals to be raced
    * seed: Seed of the racing object

    """
    # We set the defaults or the kwargs
    arg_list = []
    if(pop is not None):
        arg_list.append(pop)
    arg_list.append(seed)
    self._orig_init(*arg_list)

race_pop._orig_init = race_pop.__init__
race_pop.__init__ = _race_pop_ctor

# enum
_util.race_pop.termination_condition = _util._termination_condition


def _race_pop_run(
        self,
        n_final,
        min_trials=0,
        max_count=500,
        delta=0.05,
        racers_idx=[],
        term_cond=race_pop.termination_condition.MAX_BUDGET,
        race_best=True,
        screen_output=False):
    """
    Start a race among the individuals

    Returns a tuple of winning indices and consumed objective function evaluation.

    USAGE: race_pop.run(n_final, min_trials=0, max_count=500, delta=0.05, racers_idx=[], term_cond=util.race_pop.termination_condition.MAX_BUDGET, race_best=True, screen_output=False):

    * n_final: The desired number of winners when race ends
    * min_trials: Each individuals be evaluated at least this number of times before being dropped
    * max_count: The allow number of function evaluations (MAX_BUDGET), or the maximum number of data points to be considered for each individual (MAX_DATA_COUNT)
    * delta: Confidence level of the statistical test
    * racers_idx: Indices of the individuals to be raced, empty means to race all individuals
    * term_cond: Can be util.race_pop.termination_condition.MAX_BUDGET or util.race_pop.termination_condition.MAX_DATA_COUNT
    * race_best: When True winners are the best, otherwise winners are the worst
    * screen_output: Log racing stats at each iteration onto the screen
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(n_final)
    arg_list.append(min_trials)
    arg_list.append(max_count)
    arg_list.append(delta)
    arg_list.append(racers_idx)
    arg_list.append(term_cond)
    arg_list.append(race_best)
    arg_list.append(screen_output)
    return self._orig_run(*arg_list)

race_pop._orig_run = race_pop.run
race_pop.run = _race_pop_run


def _race_pop_size(self):
    """
    Returns the number of individuals contained in the underlying population

    USAGE: race_pop.size()
    """
    return self._orig_size()
race_pop._orig_size = race_pop.size
race_pop.size = _race_pop_size


def _race_pop_reset_cache(self):
    """
    Clears the cached fitness and constraint vectors.
    """
    return self._orig_reset_cache
race_pop._orig_reset_cache = race_pop.reset_cache
race_pop.reset_cache = _race_pop_reset_cache


def _race_pop_register_pop(self, pop):
    """Load a population into the race environment.

    This step is required before the calling to run(), if during construction
    no population was supplied.

    * pop Population to be registered. Racing will operate over this population.
    """
    return self._orig_register_pop(pop)
race_pop._orig_register_pop = race_pop.register_pop
race_pop.register_pop = _race_pop_register_pop


def _race_pop_inherit_memory(self, race_pop_src):
    """Transfer compatible evaluation history from another race_pop object.

    The source race_pop object and the current race_pop object must have the
    same seed. Upon calling, the current race_pop object will inherit
    evaluation history of individuals who also happen to reside in source.

    USAGE: race_pop.inherit_memory(race_pop_src)

    * race_pop_src: The source *race_pop* object from which compatible evaluation history will be transferred to current object
    """
    return self._orig_inherit_memory(race_pop_src)
race_pop._orig_inherit_memory = race_pop.inherit_memory
race_pop.inherit_memory = _race_pop_inherit_memory


def _race_pop_get_mean_fitness(self, ind_list=[]):
    """
    Returns the average fitness value of the individuals based on the evaluation history

    * ind_list: The indices of the individuals whose mean fitness vectors are to be extracted. If this is empty, mean data of all the individuals will be returned.
    """
    return self._orig_get_mean_fitness(ind_list)
race_pop._orig_get_mean_fitness = race_pop.get_mean_fitness
race_pop.get_mean_fitness = _race_pop_get_mean_fitness


def _race_pop_set_seed(self, seed):
    """
    Reset the seed for racing.

    * seed: The new seed to be set. This automatically clears the evaluation cache.
    """
    return self._orig_set_seed(seed)
race_pop._orig_set_seed = race_pop.set_seed
race_pop.set_seed = _race_pop_set_seed


def _race_algo_ctor(self, algo_list, probs, pop_size=100, seed=0):
    """
    Construct the racing object responsible for racing algorithms

    * algo_list: The algorithms to be raced
    * probs: Can be a single PyGMO problem or a list of them
    * pop_size: All the algorithms will be evolving internally some random population of this size
    * seed: Seed of the race
    """
    # We set the defaults or the kwargs
    arg_list = []

    #algo_vec = vector_of_algorithm_base_ptr()
    #algo_vec.extend(algo_list)
    arg_list.append(algo_list)

    try:
        l = len(probs)
        #prob_vec = vector_of_problem_base_ptr()
        #prob_vec.extend(probs)
        arg_list.append(probs)
    except TypeError:
        arg_list.append(probs)

    arg_list.append(pop_size)
    arg_list.append(seed)

    self._orig_init(*arg_list)

race_algo._orig_init = race_algo.__init__
race_algo.__init__ = _race_algo_ctor


def _race_algo_run(
        self,
        n_final,
        min_trials=0,
        max_count=500,
        delta=0.05,
        racers_idx=[],
        race_best=True,
        screen_output=False):
    """
    Start a race among several algorithms

    Returns a tuple of winning indices and the total number of evolve() made.

    USAGE: race_algo.run(n_final, min_trials=0, max_count=500, delta=0.05, racers_idx=[], race_best=True, screen_output=False):

    * n_final: The desired number of winners when race ends
    * min_trials: Each algorithms be evaluated at least this number of times before being dropped
    * max_count: The allow number of algorithm performance evaluation (i.e. number of calls to evolve)
    * delta: Confidence level of the statistical test
    * racers_idx: Indices of the algorithms to be raced, empty means to race all algorithms
    * race_best: When True winners are the best, otherwise winners are the worst
    * screen_output: Log racing stats at each iteration onto the screen
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(n_final)
    arg_list.append(min_trials)
    arg_list.append(max_count)
    arg_list.append(delta)
    arg_list.append(racers_idx)
    arg_list.append(race_best)
    arg_list.append(screen_output)
    return self._orig_run(*arg_list)

race_algo._orig_run = race_algo.run
race_algo.run = _race_algo_run
