# -*- coding: utf-8 -*-
from PyGMO.algorithm._algorithm import *
from PyGMO.algorithm._algorithm import _base
from PyGMO.algorithm._base import base
from PyGMO.algorithm._example import py_example
from PyGMO.algorithm._cmaes import py_cmaes
from PyGMO.algorithm._scipy_algos import *

_base = _algorithm._base

# Creating the list of algorithms


def _get_algorithm_list():
    from PyGMO.algorithm import _algorithm as algorithm
    # Try importing SciPy and NumPy.
    try:
        import scipy
        import numpy
        algorithm_list = [
            algorithm.__dict__[n]
            for n in
            [n
             for n in dir(algorithm)
             if not n.startswith('_') and not n == 'base']]
    except ImportError as e:
        algorithm_list = [
            algorithm.__dict__[n]
            for n in
            [n
             for n in dir(algorithm)
             if not n.startswith('_') and not n == 'base' and not n.
             startswith('scipy')]]
    return algorithm_list


# Redefining the constructors of all algorithms to obtain good
# documentation and to allow kwargs
def _de_ctor(
        self,
        gen=100,
        f=0.8,
        cr=0.9,
        variant=2,
        ftol=1e-6,
        xtol=1e-6,
        screen_output=False):
    """
    Constructs a Differential Evolution algorithm:

    USAGE: algorithm.de(gen=1, f=0.5, cr=0.9, variant=2, ftol=1e-6, xtol=1e-6, screen_output = False)

    * gen: number of generations
    * f: weighting factor in [0,1] (if -1 self-adptation is used)
    * cr: crossover in [0,1] (if -1 self-adptation is used)
    * variant: algoritmic variant to use (one of [1 .. 10])
            1. DE/best/1/exp
            2. DE/rand/1/exp
            3. DE/rand-to-best/1/exp
            4. DE/best/2/exp
            5. DE/rand/2/exp
            6. DE/best/1/bin
            7. DE/rand/1/bin
            8. DE/rand-to-best/1/bin
            9. DE/best/2/bin
            10. DE/rand/2/bin
    * ftol stop criteria on f
    * xtol stop criteria on x
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(gen)
    arg_list.append(f)
    arg_list.append(cr)
    arg_list.append(variant)
    arg_list.append(ftol)
    arg_list.append(xtol)
    self._orig_init(*arg_list)
    self.screen_output = screen_output
de._orig_init = de.__init__
de.__init__ = _de_ctor


def _jde_ctor(
        self,
        gen=100,
        variant=2,
        variant_adptv=1,
        ftol=1e-6,
        xtol=1e-6,
        memory=False,
        screen_output=False):
    """
    Constructs a jDE algorithm (self-adaptive DE)

    REF: "Self-adaptive differential evolution algorithm in constrained real-parameter optimization"
    J Brest, V Zumer, MS Maucec - Evolutionary Computation, 2006.
    http://dsp.szu.edu.cn/DSP2006/research/publication/yan/WebEdit/UploadFile/Self-adaptive%20Differential%20Evolution%20Algorithm%20for%20Constrained%20Real-Parameter%20Optimization.pdf

    USAGE: algorithm.jde(gen=100, variant=2, variant_adptv=1, ftol=1e-6, xtol=1e-6, memory = False, screen_output = False)

    * gen: number of generations
    * variant: algoritmic variant to use (one of [1 .. 18])
            1. best/1/exp				2. rand/1/exp
            3. rand-to-best/1/exp		4. best/2/exp
            5. rand/2/exp				6. best/1/bin
            7. rand/1/bin				8. rand-to-best/1/bin
            9. best/2/bin				10. rand/2/bin
            11. best/3/exp				12. best/3/bin
            13. rand/3/exp				14. rand/3/bin
            15. rand-to-current/2/exp		16. rand-to-current/2/bin
            17. rand-to-best-and-current/2/exp	18. rand-to-best-and-current/2/bin
    * variant_adptv: adaptive scheme to use (one of [1..2])
            1. random param mutation		2. param mutation follows rand/3 scheme
    * ftol: stop criteria on f
    * xtol: stop criteria on x
    * memory: if True the algorithm internal state is saved and used for the next call
    * screen_output: activates screen output of the algorithm (do not use in archipealgo, otherwise the screen will be flooded with
    *				 different island outputs)
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(gen)
    arg_list.append(variant)
    arg_list.append(variant_adptv)
    arg_list.append(ftol)
    arg_list.append(xtol)
    arg_list.append(memory)
    self._orig_init(*arg_list)
    self.screen_output = screen_output
jde._orig_init = jde.__init__
jde.__init__ = _jde_ctor


def _de_1220_ctor(
        self,
        gen=100,
        variant_adptv=1,
        allowed_variants=[
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10],
        memory=False,
        ftol=1e-6,
        xtol=1e-6,
        screen_output=False):
    """
    Constructs a Differential Evolution algorithm (our own brew). Self adaptation on F, CR and mutation variant.:

    USAGE: algorithm.de_1220(gen=100, variant_adptv=1, allowed_variants = [i for i in range(1,19)], memory = False, ftol=1e-6, xtol=1e-6, screen_output = False)

    * gen: number of generations
    * variant_adptv: adaptiv scheme to use (one of [1..2])
            1. random param mutation		2. param mutation follows relative DE scheme
    * allowed_variants : a list of the algoritmic variants to mix and self-adapt. Allowed variants are ...
            1. best/1/exp				2. rand/1/exp
            3. rand-to-best/1/exp			4. best/2/exp
            5. rand/2/exp				6. best/1/bin
            7. rand/1/bin				8. rand-to-best/1/bin
            9. best/2/bin				10. rand/2/bin
            11. best/3/exp				12. best/3/bin
            13. rand/3/exp				14. rand/3/bin
            15. rand-to-current/2/exp		16. rand-to-current/2/bin
            17. rand-to-best-and-current/2/exp	18. rand-to-best-and-current/2/bin
    * ftol: stop criteria on f
    * xtol: stop criteria on x
    * memory: if True the algorithm internal state is saved and used for the next call
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(gen)
    arg_list.append(variant_adptv)
    arg_list.append(allowed_variants)
    arg_list.append(memory)
    arg_list.append(ftol)
    arg_list.append(xtol)
    self._orig_init(*arg_list)
    self.screen_output = screen_output
de_1220._orig_init = de_1220.__init__
de_1220.__init__ = _de_1220_ctor


def _mde_pbx_ctor(
        self,
        gen=100,
        qperc=0.15,
        nexp=1.5,
        ftol=1e-6,
        xtol=1e-6,
        screen_output=False):
    """
    Constructs a mde_pbx algorithm (self-adaptive DE)

    REF: "An Adaptive Differential Evolution Algorithm With Novel Mutation and Crossover
    Strategies for Global Numerical Optimization" - IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS?PART B: CYBERNETICS, VOL. 42, NO. 2, APRIL 20


    USAGE: algorithm.mde_pbx(gen=100, qperc=0.15, nexp=1.5, ftol=1e-6, xtol=1e-6, screen_output = False)

    * gen: number of generations
    * qperc: percentage of population to choose the best vector
    * nexp: exponent for the powermean
    * ftol: stop criteria on f
    * xtol: stop criteria on x
    * screen_output: activates screen output of the algorithm (do not use in archipealgo, otherwise the screen will be flooded with
    * 		 different island outputs)
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(gen)
    arg_list.append(qperc)
    arg_list.append(nexp)
    arg_list.append(ftol)
    arg_list.append(xtol)
    self._orig_init(*arg_list)
    self.screen_output = screen_output
mde_pbx._orig_init = mde_pbx.__init__
mde_pbx.__init__ = _mde_pbx_ctor


def _pso_ctor(
        self,
        gen=1,
        omega=0.7298,
        eta1=2.05,
        eta2=2.05,
        vcoeff=0.5,
        variant=5,
        neighb_type=2,
        neighb_param=4):
    """
    Constructs a Particle Swarm Optimization (steady-state). The position update is applied
    immediately after the velocity update

    REF (for variants 5-6): http://cswww.essex.ac.uk/staff/rpoli/papers/PoliKennedyBlackwellSI2007.pdf

    REF (for variants 1-4): Kennedy, J.; Eberhart, R. (1995). "Particle Swarm Optimization". Proceedings of IEEE International Conference on Neural Networks. IV. pp. 1942?1948.

    USAGE: algorithm.pso(gen=1, omega = 0.7298, eta1 = 2.05, eta2 = 2.05, vcoeff = 0.5, variant = 5, neighb_type = 2, neighb_param = 4)

    * gen: number of generations
    * omega: constriction factor (or particle inertia weight) in [0,1]
    * eta1: Cognitive component in [0,4]
    * eta2: Social component in [0,4]
    * vcoeff: Maximum velocity coefficient (w.r.t. the box-bounds width) in [0,1]
    * variant: algoritmic variant to use (one of  [1 .. 6])
            1. PSO canonical (with inertia weight)
            2. PSO canonical (with inertia weight
                    and equal random weights of social and cognitive components)
            3. PSO variant (with inertia weight
                    same random number for all components.)
            4. PSO variant (with inertia weight
                    same random number for all components
                    and equal weights of social and cognitive components)
            5. PSO canonical (with constriction factor)
            6. Fully Informed Particle Swarm (FIPS)
    * neighb_type: defines the particle neighbourhood (used for the social component)
            1. gbest neighbourhood topology (fully connected)
            2. lbest neighbourhood topology (ring)
            3. Von-Neumann neighbourhood topology (square lattice)
            4. Randomly-varying neighbourhood topology
    * neighb_param: if the lbest topology is selected, it represents each particle's indegree
            (also outdegree) in the swarm topology. Particles have neighbours up
            to a radius of k = neighb_param / 2 in the ring. If the Randomly-varying neighbourhood topology
            is selected, neighb_param represents each particle's maximum outdegree in the swarm topology.
            The minimum outdegree is 1 (the particle always connects back to itself).
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(gen)
    arg_list.append(omega)
    arg_list.append(eta1)
    arg_list.append(eta2)
    arg_list.append(vcoeff)
    arg_list.append(variant)
    arg_list.append(neighb_type)
    arg_list.append(neighb_param)
    self._orig_init(*arg_list)
pso._orig_init = pso.__init__
pso.__init__ = _pso_ctor


def _pso_gen_ctor(
        self,
        gen=1,
        omega=0.7298,
        eta1=2.05,
        eta2=2.05,
        vcoeff=0.5,
        variant=5,
        neighb_type=2,
        neighb_param=4):
    """
    Constructs a Particle Swarm Optimization (generational). The position update is applied
    only at the end of an entire loop over the population (swarm). Use this version for stochastic problems.

    USAGE: algorithm.pso_gen(gen=1, omega = 0.7298, eta1 = 2.05, eta2 = 2.05, vcoeff = 0.5, variant = 5, neighb_type = 2, neighb_param = 4)

    * gen: number of generations
    * omega: constriction factor (or particle inertia weight) in [0,1]
    * eta1: Cognitive component in [0,4]
    * eta2: Social component in [0,4]
    * vcoeff: Maximum velocity coefficient (w.r.t. the box-bounds width) in [0,1]
    * variant: algoritmic variant to use (one of  [1 .. 6])
            1. PSO canonical (with inertia weight)
            2. PSO canonical (with inertia weight
                    and equal random weights of social and cognitive components)
            3. PSO variant (with inertia weight
                    same random number for all components.)
            4. PSO variant (with inertia weight
                    same random number for all components
                    and equal weights of social and cognitive components)
            5. PSO canonical (with constriction factor)
            6. Fully Informed Particle Swarm (FIPS)
    * neighb_type: defines the particle neighbourhood (used for the social component)
            1. gbest neighbourhood topology (fully connected)
            2. lbest neighbourhood topology (ring)
            3. Von-Neumann neighbourhood topology (square lattice)
            4. Randomly-varying neighbourhood topology
    * neighb_param: if the lbest topology is selected, it represents each particle's indegree
            (also outdegree) in the swarm topology. Particles have neighbours up
            to a radius of k = neighb_param / 2 in the ring. If the Randomly-varying neighbourhood topology
            is selected, neighb_param represents each particle's maximum outdegree in the swarm topology.
            The minimum outdegree is 1 (the particle always connects back to itself).
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(gen)
    arg_list.append(omega)
    arg_list.append(eta1)
    arg_list.append(eta2)
    arg_list.append(vcoeff)
    arg_list.append(variant)
    arg_list.append(neighb_type)
    arg_list.append(neighb_param)
    self._orig_init(*arg_list)
pso_gen._orig_init = pso_gen.__init__
pso_gen.__init__ = _pso_gen_ctor


def _pso_gen_racing_ctor(
        self,
        gen=1,
        omega=0.7298,
        eta1=2.05,
        eta2=2.05,
        vcoeff=0.5,
        variant=5,
        neighb_type=2,
        neighb_param=4,
        nr_eval_per_x=5,
        max_fevals=10000000):
    """
    Constructs a Particle Swarm Optimization (generational). The position update is applied
    only at the end of an entire loop over the population (swarm). Use this version for stochastic problems.

    USAGE: algorithm.pso_gen(gen=1, omega = 0.7298, eta1 = 2.05, eta2 = 2.05, vcoeff = 0.5, variant = 5, neighb_type = 2, neighb_param = 4, nr_eval_per_x = 5, max_fevals = 10000000)

    * gen: number of generations
    * omega: constriction factor (or particle inertia weight) in [0,1]
    * eta1: Cognitive component in [0,4]
    * eta2: Social component in [0,4]
    * vcoeff: Maximum velocity coefficient (w.r.t. the box-bounds width) in [0,1]
    * variant: algoritmic variant to use (one of  [1 .. 6])
            1. PSO canonical (with inertia weight)
            2. PSO canonical (with inertia weight
                    and equal random weights of social and cognitive components)
            3. PSO variant (with inertia weight
                    same random number for all components.)
            4. PSO variant (with inertia weight
                    same random number for all components
                    and equal weights of social and cognitive components)
            5. PSO canonical (with constriction factor)
            6. Fully Informed Particle Swarm (FIPS)
    * neighb_type: defines the particle neighbourhood (used for the social component)
            1. gbest neighbourhood topology (fully connected)
            2. lbest neighbourhood topology (ring)
            3. Von-Neumann neighbourhood topology (square lattice)
            4. Randomly-varying neighbourhood topology
    * neighb_param: if the lbest topology is selected, it represents each particle's indegree
            (also outdegree) in the swarm topology. Particles have neighbours up
            to a radius of k = neighb_param / 2 in the ring. If the Randomly-varying neighbourhood topology
            is selected, neighb_param represents each particle's maximum outdegree in the swarm topology.
            The minimum outdegree is 1 (the particle always connects back to itself).
* nr_eval_per_x: Specify the expected budget to be allocated during racing
    * max_fevals: When specified other than -1, this serve as another termination condition -- maximium number of objective function evaluations
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(gen)
    arg_list.append(omega)
    arg_list.append(eta1)
    arg_list.append(eta2)
    arg_list.append(vcoeff)
    arg_list.append(variant)
    arg_list.append(neighb_type)
    arg_list.append(neighb_param)
    arg_list.append(nr_eval_per_x)
    if max_fevals > 0:
        arg_list.append(max_fevals)
    self._orig_init(*arg_list)
pso_gen_racing._orig_init = pso_gen_racing.__init__
pso_gen_racing.__init__ = _pso_gen_racing_ctor

_algorithm.sga.crossover = _algorithm._sga_crossover_type
_algorithm.sga.selection = _algorithm._sga_selection_type
_algorithm.sga.mutation = _algorithm._sga_mutation_type


def _sga_ctor(
        self,
        gen=1,
        cr=.95,
        m=.02,
        elitism=1,
        mutation=sga.mutation.GAUSSIAN,
        width=0.1,
        selection=sga.selection.ROULETTE,
        crossover=sga.crossover.EXPONENTIAL):
    """
    Constructs a Simple Genetic Algorithm (generational)

    USAGE: algorithm.sga(self, gen=1, cr=.95, m=.02, elitism=1, mutation=sga.mutation.GAUSSIAN, width = 0.1, selection=sga.selection.ROULETTE, crossover=sga.crossover.EXPONENTIAL)

    * gen: number of generations
    * cr: crossover factor in [0,1]
    * m: mutation probability (for each component) [0,1]
    * elitism: number of generation after which the best is reinserted
    * mutation: mutation type (one of [RANDOM, GAUSSIAN])
    * width: the mutation width (in case of a GAUSSIAN bell
            this is the std normalized with the width)
    * selection: selection startegy (one of [ROULETTE, BEST20])
    * crossover: crossover strategy (one of [BINOMIAL, EXPONENTIAL])
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(gen)
    arg_list.append(cr)
    arg_list.append(m)
    arg_list.append(elitism)
    arg_list.append(mutation)
    arg_list.append(width)
    arg_list.append(selection)
    arg_list.append(crossover)
    self._orig_init(*arg_list)
sga._orig_init = sga.__init__
sga.__init__ = _sga_ctor

_algorithm.vega.crossover = _algorithm._vega_crossover_type
_algorithm.vega.mutation = _algorithm._vega_mutation_type


def _vega_ctor(
        self,
        gen=1,
        cr=.95,
        m=.02,
        elitism=1,
        mutation=vega.mutation.GAUSSIAN,
        width=0.1,
        crossover=vega.crossover.EXPONENTIAL):
    """
    Constructs a Vector evaluated genetic algorithm

    USAGE: algorithm.vega(self, gen=1, cr=.95, m=.02, elitism=1, mutation=vega.mutation.GAUSSIAN, width = 0.1, crossover=vega.crossover.EXPONENTIAL)

    * gen: number of generations
    * cr: crossover factor in [0,1]
    * m: mutation probability (for each component) [0,1]
    * elitism: number of generation after which the best is reinserted
    * mutation: mutation type (one of [RANDOM, GAUSSIAN])
    * width: the mutation width (in case of a GAUSSIAN bell
            this is the std normalized with the width)
    * crossover: crossover strategy (one of [BINOMIAL, EXPONENTIAL])
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(gen)
    arg_list.append(cr)
    arg_list.append(m)
    arg_list.append(elitism)
    arg_list.append(mutation)
    arg_list.append(width)
    arg_list.append(crossover)
    self._orig_init(*arg_list)
vega._orig_init = vega.__init__
vega.__init__ = _vega_ctor

_algorithm.sga_gray.crossover = _algorithm._gray_crossover_type
_algorithm.sga_gray.selection = _algorithm._gray_selection_type
_algorithm.sga_gray.mutation = _algorithm._gray_mutation_type


def _sga_gray_ctor(
        self,
        gen=1,
        cr=.95,
        m=.02,
        elitism=1,
        mutation=sga_gray.mutation.UNIFORM,
        selection=sga_gray.selection.ROULETTE,
        crossover=sga_gray.crossover.SINGLE_POINT):
    """
    Constructs a Simple Genetic Algorithm with gray binary encoding (generational)

    USAGE: algorithm.sga_gray(self, gen=1, cr=.95, m=.02, elitism=1, mutation=sga.mutation.UNIFORM, selection=sga.selection.ROULETTE, crossover=sga.crossover.SINGLE_POINT)

    * gen: Number of generations to evolve.
    * cr: crossover factor in [0,1]
    * m: mutation probability (of each encoded bit) [0,1]
    * elitism: number of generation after which the best is reinserted
    * mut: mutation type (one of [UNIFORM])
    * sel: selection strategy (one of [ROULETTE, BEST20])
    * cro: crossover strategy (one of [SINGLE_POINT])
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(gen)
    arg_list.append(cr)
    arg_list.append(m)
    arg_list.append(elitism)
    arg_list.append(mutation)
    arg_list.append(selection)
    arg_list.append(crossover)
    self._orig_init(*arg_list)
sga_gray._orig_init = sga_gray.__init__
sga_gray.__init__ = _sga_gray_ctor


def _nsga_II_ctor(self, gen=100, cr=0.95, eta_c=10, m=0.01, eta_m=10):
    """
    Constructs a Non-dominated Sorting Genetic Algorithm (NSGA_II)

    USAGE: algorithm.nsga_II(self, gen=100, cr = 0.95, eta_c = 10, m = 0.01, eta_m = 10)

    * gen: number of generations
    * cr: crossover factor [0,1[
    * eta_c: Distribution index for crossover
    * m: mutation probability [0,1]
    * eta_m: Distribution index for mutation
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(gen)
    arg_list.append(cr)
    arg_list.append(eta_c)
    arg_list.append(m)
    arg_list.append(eta_m)
    self._orig_init(*arg_list)
nsga_II._orig_init = nsga_II.__init__
nsga_II.__init__ = _nsga_II_ctor


def _sms_emoa_ctor(
        self,
        hv_algorithm=None,
        gen=100,
        sel_m=2,
        cr=0.95,
        eta_c=10,
        m=0.01,
        eta_m=10):
    """
    Constructs a S-Metric Selection Evolutionary Multiobjective Optimiser Algorithm (SMS-EMOA)

    USAGE: algorithm.sms_emoa(self, gen=100, sel_m = 2, cr = 0.95, eta_c = 10, m = 0.01, eta_m = 10)

    * hv_algorithm: hypervolume algorithm object used for the computation of the hypervolume. By default its chosen dynamically
    * gen: number of generations
    * sel_m: selection method for points in dominated fronts. 1 - always use least contributor, 2 - use domination count for fronts > 1
    * cr: crossover factor [0,1]
    * eta_c: Distribution index for crossover
    * m: mutation probability [0,1]
    * eta_m: Distribution index for mutation
    """
    # We set the defaults or the kwargs
    arg_list = []

    from PyGMO.util import _HypervolumeValidation
    if hv_algorithm:
        hv_algorithm = _HypervolumeValidation.validate_hv_algorithm(
            hv_algorithm)
        arg_list.append(hv_algorithm)
    arg_list.append(gen)
    arg_list.append(sel_m)
    arg_list.append(cr)
    arg_list.append(eta_c)
    arg_list.append(m)
    arg_list.append(eta_m)
    self._orig_init(*arg_list)
sms_emoa._orig_init = sms_emoa.__init__
sms_emoa.__init__ = _sms_emoa_ctor


def _pade_ctor(
        self,
        gen=10,
        decomposition='tchebycheff',
        weights='grid',
        solver=None,
        threads=8,
        T=8,
        z=[]):
    """
    Constructs a Parallel Decomposition Algorithm (PaDe).

    For each element of the population a different single objective problem is generated using a decomposition method.
    Those single-objective problems are thus solved in an island model.
    At the end of the evolution the population is set as the best individual in each single-objective island.
    This algorithm, original with PaGMO, builds upon the MOEA/D framework

    USAGE: algorithm.pade(self, gen=10, , decomposition = 'tchebycheff', weights = 'grid', solver = None, threads = 8, T = 8, z = [])

    * gen: number of generations
    * threads: the maximum number of single-objective problems to solve at the same time
    * solver: the algorithm to use to solve the single-objective problems
    * T: the size of the population on each subproblem (must be an even number)
    * decomposition = the decomposition method to use, on of  ('weighted', 'tchebycheff' or 'bi')
    * weights: weight generation method, one of ('grid', 'low_discrepancy', 'random')
    * z: the reference point (used with Tchebycheff and BI decomposition methods)
    """
    # We set the defaults or the kwargs
    from PyGMO.problem._problem_meta import _decomposition_method
    DECOMPOSITION_TYPE = {
        'weighted': _decomposition_method.WEIGHTED,
        'tchebycheff': _decomposition_method.TCHEBYCHEFF,
        'bi': _decomposition_method.BI,
    }

    WEIGHT_GENERATION_TYPE = {
        'low_discrepancy': _algorithm._weight_generation.LOW_DISCREPANCY,
        'grid': _algorithm._weight_generation.GRID,
        'random': _algorithm._weight_generation.RANDOM,
    }

    arg_list = []
    arg_list.append(gen)
    arg_list.append(threads)
    arg_list.append(DECOMPOSITION_TYPE[decomposition.lower()])
    if solver is None:
        solver = jde(100)
    arg_list.append(solver)
    arg_list.append(T)
    arg_list.append(WEIGHT_GENERATION_TYPE[weights.lower()])
    arg_list.append(z)
    self._orig_init(*arg_list)
pade._orig_init = pade.__init__
pade.__init__ = _pade_ctor


def _nspso_ctor(
        self,
        gen=100,
        minW=0.4,
        maxW=1.0,
        C1=2.0,
        C2=2.0,
        CHI=1.0,
        v_coeff=0.5,
        leader_selection_range=5,
        diversity_mechanism='crowding distance'):
    """
    Constructs a Multi Objective PSO

    USAGE: algorithm.nspso(self, gen=10, minW = 0.4, maxW = 1.0, C1 = 2.0, C2 = 2.0,
            CHI = 1.0, v_coeff = 0.5, leader_selection = 5, diversity_mechanism = 'crowding_distance'):

    * gen: number of generations
    * minW: minimum particles' inertia weight (the inertia weight is decreased troughout the run between maxW and minW)
    * maxW: maximum particles' inertia weight (the inertia weight is decreased troughout the run between maxW and minW)
    * C1: magnitude of the force, applied to the particle's velocity, in the direction of its previous best position
    * C2: magnitude of the force, applied to the particle's velocity, in the direction of its global best (leader)
    * CHI: velocity scaling factor
    * v_coeff: velocity coefficient (determining the maximum allowed particle velocity)
    * leader_selection_range the leader of each particle is selected among the best leader_selection_range% individuals
    * diversity_mechanism one of 'crowding distance', 'niche count' or 'maxmin'. Defines the diversity mechanism used
    """
    # We set the defaults or the kwargs
    DIVERSITY_MECHANISM = {
        'crowding distance': _algorithm._diversity_mechanism.CROWDING_DISTANCE,
        'niche count': _algorithm._diversity_mechanism.NICHE_COUNT,
        'maxmin': _algorithm._diversity_mechanism.MAXMIN,
    }

    arg_list = []
    arg_list.append(gen)
    arg_list.append(minW)
    arg_list.append(maxW)
    arg_list.append(C1)
    arg_list.append(C2)
    arg_list.append(CHI)
    arg_list.append(v_coeff)
    arg_list.append(leader_selection_range)
    arg_list.append(DIVERSITY_MECHANISM[diversity_mechanism])
    self._orig_init(*arg_list)

nspso._orig_init = nspso.__init__
nspso.__init__ = _nspso_ctor


def _spea2_ctor(
        self,
        gen=100,
        cr=0.95,
        eta_c=10,
        m=0.01,
        eta_m=50,
        archive_size=0):
    """
    Constructs a Strenght Pareto Evolutionary Algorithm 2

    USAGE: algorithm.spea2(gen=100, cr = 0.95, eta_c = 10, m = 0.01, eta_m = 50, archive_size = -1)

    * gen: Number of generations to evolve.
    * cr: Crossover probability
    * eta_c: Distribution index for crossover
    * m: Mutation probability
    * eta_m: Distribution index for mutation
    * archive_size: the size of the non_dominated archive. If archive_size=0 then the archive size is set equal to the population size. The population returned after evolve has a size equal to archive_size
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(gen)
    arg_list.append(cr)
    arg_list.append(eta_c)
    arg_list.append(m)
    arg_list.append(eta_m)
    arg_list.append(archive_size)
    self._orig_init(*arg_list)

spea2._orig_init = spea2.__init__
spea2.__init__ = _spea2_ctor

from PyGMO.problem import decompose


def _moead_ctor(
        self,
        gen=100,
        weights='grid',
        T=20,
        realb=0.9,
        limit=2,
        cr=1.0,
        f=0.5,
        eta_m=20,
        diversity=True):
    """
    Multi Objective Evolutionary Algorithm based on Decomposition and Differential Evolution (MOEA/D - DE)

    REF Zhang, Qingfu, and Hui Li. "MOEA/D: A multiobjective evolutionary algorithm based on decomposition." Evolutionary Computation, IEEE Transactions on 11.6 (2007): 712-731.
    REF Li, Hui, and Qingfu Zhang. "Multiobjective optimization problems with complicated Pareto sets, MOEA/D and NSGA-II." Evolutionary Computation, IEEE Transactions on 13.2 (2009): 284-302.

    USAGE: algorithm.spea2(gen=100, weights = 'grid', T = 20, realb = 0.9, limit = 2, cr = 1.0, f = 0.5, eta_m=20, diversity=True)

    * gen: Number of generations to evolve.
    * weights: weight generation method, one of ('grid', 'low_discrepancy', 'random')
    * T: Size of the neighbourhood
    * realb Chance that the neighbourhood is T rather than the whole population (only if diversity is True)
    * limit Maximum number of copies reinserted in the population  (only if diversity is True)
    * cr Crossover parameter in the Differential Evolution operator
    * f f parameter in the Differential Evolution operator
    * eta_m Distribution index for the polynomial mutation
    * diversity when true activates the two diversity preservation mechanism described in Li, Hui, and Qingfu Zhang paper
    """
    def weight_generation_type(x):
        return {
            'low_discrepancy': _algorithm._weight_generation_moead.LOW_DISCREPANCY,
            'grid': _algorithm._weight_generation_moead.GRID,
            'random': _algorithm._weight_generation_moead.RANDOM,
        }[x]

    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(gen)
    arg_list.append(weight_generation_type(weights))
    arg_list.append(T)
    arg_list.append(realb)
    arg_list.append(limit)
    arg_list.append(cr)
    arg_list.append(f)
    arg_list.append(eta_m)
    arg_list.append(diversity)
    self._orig_init(*arg_list)
moead._orig_init = moead.__init__
moead.__init__ = _moead_ctor


def _sa_corana_ctor(
        self,
        iter=10000,
        Ts=10,
        Tf=.1,
        steps=1,
        bin_size=20,
        range=1):
    """
    Constructs Corana's Simulated Annealing

    USAGE: algorithm.sa_corana(iter = 10000, Ts = 10, Tf = .1, steps = 1, bin_size = 20, range = 1)

    NOTE: as this version of simulated annealing loops through the chromosome, the iter number needs to be selected
    large enough to allow the temperature schedule to actuallt make sense. For example if your problem has D dimensions
    then in order to have at least N temperature adjustments (from Ts to Tf) one should select iter = D * N * steps * bin_size.

    * iter: number of total iterations
    * Ts: starting temperature
    * Tf: final temperature ( > Ts)
    * steps: number of steps adjustments
    * bin_size: size of the bin used to evaluate the step adjustment
    * range: initial size of the neighbourhood (in [0,1])
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(iter)
    arg_list.append(Ts)
    arg_list.append(Tf)
    arg_list.append(steps)
    arg_list.append(bin_size)
    arg_list.append(range)
    self._orig_init(*arg_list)
sa_corana._orig_init = sa_corana.__init__
sa_corana.__init__ = _sa_corana_ctor


def _bee_colony_ctor(self, gen=100, limit=20):
    """
    Constructs an Artificial Bee Colony Algorithm

    USAGE: algorithm.bee_colony(gen = 100, limit = 20)

    * gen: number of 'generations' (each generation 2*NP function evaluations
            are made where NP is the population size)
    * limit: number of tries after which a source of food is dropped if not improved
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(gen)
    arg_list.append(limit)
    self._orig_init(*arg_list)
bee_colony._orig_init = bee_colony.__init__
bee_colony.__init__ = _bee_colony_ctor


def _sea_ctor(self, gen=100, limit=20):
    """
    Constructs a simple (N+1)-EA: A Simple Evolutionary Algorithm

    USAGE: algorithm.ea(gen = 1)
    SEE : Oliveto, Pietro S., Jun He, and Xin Yao.
    "Time complexity of evolutionary algorithms for combinatorial optimization: A decade of results."
    International Journal of Automation and Computing 4.3 (2007): 281-293.

    * gen: number of 'generations' (each generation is one function evaluation)
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(gen)
    self._orig_init(*arg_list)
sea._orig_init = sea.__init__
sea.__init__ = _sea_ctor


def _ms_ctor(self, algorithm=None, iter=1):
    """
    Constructs a Multistart Algorithm

    USAGE: algorithm.ms(algorithm = algorithm.de(), iter = 1)

    NOTE: starting from pop1, at each iteration a random pop2 is evolved
    with the selected algorithm and its final best replaces the worst of pop1

    * algorithm: PyGMO algorithm to be multistarted
    * iter: number of multistarts

    """
    # We set the defaults or the kwargs
    arg_list = []
    if algorithm is None:
        algorithm = _algorithm.jde()
    arg_list.append(algorithm)
    arg_list.append(iter)
    self._orig_init(*arg_list)
ms._orig_init = ms.__init__
ms.__init__ = _ms_ctor


def _cs_ctor(
        self,
        max_eval=1,
        stop_range=0.01,
        start_range=0.1,
        reduction_coeff=0.5):
    """
    Constructs a Compass Search Algorithm

    USAGE: algorithm.cs(max_eval = 1, stop_range = 0.01, start_range = 0.1, reduction_coeff = 0.5);


    * max_eval: maximum number of function evaluations
    * stop_range: when the range is reduced to a value smaller than stop_range cs stops
    * start_range: starting range (non-dimensional wrt ub-lb)
    * reduction_coeff: the range is multiplied by reduction_coeff whenever no improvment is made
                       across one chromosome
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(max_eval)
    arg_list.append(stop_range)
    arg_list.append(start_range)
    arg_list.append(reduction_coeff)
    self._orig_init(*arg_list)
cs._orig_init = cs.__init__
cs.__init__ = _cs_ctor


def _mbh_ctor(self, algorithm=None, stop=5, perturb=5e-2, screen_output=False):
    """
    Constructs a Monotonic Basin Hopping Algorithm (generalized to accept any algorithm)

    USAGE: algorithm.mbh(algorithm = algorithm.cs(), stop = 5, perturb = 5e-2);

    NOTE: Starting from pop, algorithm is applied to the perturbed pop returning pop2. If pop2 is better than
    pop then pop=pop2 and a counter is reset to zero. If pop2 is not better the counter is incremented. If
    the counter is larger than stop, optimization is terminated

    * algorithm: 'local' optimiser
    * stop: number of no improvements before halting the optimization
    * perturb: non-dimentional perturbation width (can be a list, in which case
            it has to have the same dimension of the problem mbh will be applied to)
    * screen_output: activates screen output of the algorithm (do not use in archipealgo, otherwise the screen will be flooded with
    * 		 different island outputs)
    """
    # We set the defaults or the kwargs
    arg_list = []
    if algorithm is None:
        algorithm = _algorithm.cs()
    arg_list.append(algorithm)
    arg_list.append(stop)
    arg_list.append(perturb)
    self._orig_init(*arg_list)
    self.screen_output = screen_output
mbh._orig_init = mbh.__init__
mbh.__init__ = _mbh_ctor


def _cstrs_self_adaptive_ctor(
        self,
        algorithm=None,
        max_iter=100,
        f_tol=1e-15,
        x_tol=1e-15):
    """
    Constructs a Self-Adaptive Fitness constraints handling Meta Algorithm.

    The key idea of this constraint handling technique is to represent the
    constraint violation by a single infeasibility measure, and to adapt
    dynamically the penalization of infeasible solutions.

    USAGE: algorithm.self_adaptive(algorithm = algorithm.jde(), max_iter  = 100, f_tol = 1e-15, x_tol = 1e-15);

    * algorithm: original optimizer
    * max_iter: stop-criteria (number of iterations)
    * ftol: 1e-15 by default. The stopping criteria on the x tolerance.
    * xtol: 1e-15 by default. The stopping criteria on the f tolerance.
    """
    # We set the defaults or the kwargs
    arg_list = []
    if algorithm is None:
        algorithm = _algorithm.jde()
    arg_list.append(algorithm)
    arg_list.append(max_iter)
    arg_list.append(f_tol)
    arg_list.append(x_tol)
    self._orig_init(*arg_list)
cstrs_self_adaptive._orig_init = cstrs_self_adaptive.__init__
cstrs_self_adaptive.__init__ = _cstrs_self_adaptive_ctor

# Renaming and placing the enums
_algorithm.cstrs_co_evolution.method = _algorithm._co_evo_method_type


def _cstrs_co_evolution_ctor(
        self,
        original_algo=None,
        original_algo_penalties=None,
        pop_penalties_size=30,
        gen=20,
        method=cstrs_co_evolution.method.SIMPLE,
        pen_lower_bound=0.,
        pen_upper_bound=100000.,
        f_tol=1e-15,
        x_tol=1e-15):
    """
    Constructs a co-evolution adaptive penalty algorithm for constrained optimization.

    USAGE: algorithm.cstrs_co_evolution(original_algo = _algorithm.jde(), original_algo_penalties = _algorithm.jde(), pop_penalties_size = 30, gen = 20, method = cstrs_co_evolution.method.SIMPLE, pen_lower_bound = 0, pen_upper_bound = 100000,f_tol = 1e-15,x_tol = 1e-15):

    * original_algo: optimizer to use as 'original' optimization method
    * original_algo_penalties: optimizer to use as 'original' optimization method for population encoding penalties coefficients
    * pop_penalties_size: size of the population encoding the penalty parameters.
    * gen: number of generations.
    * method: cstrs_co_evolution.method.SIMPLE by default, the method used for the population encoding penalties coefficients.
            Three posssibililties are available: SIMPLE,
            SPLIT_NEQ_EQ and SPLIT_CONSTRAINTS. The simple one is the original
            version of the Coello/He implementation (one penalty coefficient weights
            the sum of the constraints violation, one the number of violated constraints).
            The SPLIT_NEQ_EQ, splits the equalities and inequalities constraints in two different sets for the
            penalty weigths, containing respectively inequalities and equalities
            weigths. The SPLIT_CONSTRAINTS splits the constraints in M set of weigths
            with M the number of constraints.
    * pen_lower_bound: the lower boundary used for penalty.
    * pen_upper_bound: the upper boundary used for penalty.
    * ftol: 1e-15 by default. The stopping criteria on the x tolerance.
    * xtol: 1e-15 by default. The stopping criteria on the f tolerance.
    """
    arg_list = []
    if original_algo is None:
        original_algo = _algorithm.jde()
    if original_algo_penalties is None:
        original_algo_penalties = _algorithm.jde()
    arg_list.append(original_algo)
    arg_list.append(original_algo_penalties)
    arg_list.append(pop_penalties_size)
    arg_list.append(gen)
    arg_list.append(method)
    arg_list.append(pen_lower_bound)
    arg_list.append(pen_upper_bound)
    arg_list.append(f_tol)
    arg_list.append(x_tol)
    self._orig_init(*arg_list)
cstrs_co_evolution._orig_init = cstrs_co_evolution.__init__
cstrs_co_evolution.__init__ = _cstrs_co_evolution_ctor

# Renaming and placing the enums
_algorithm.cstrs_immune_system.select_method = _algorithm._immune_select_method_type
_algorithm.cstrs_immune_system.inject_method = _algorithm._immune_inject_method_type
_algorithm.cstrs_immune_system.distance_method = _algorithm._immune_distance_method_type


def _cstrs_immune_system_ctor(
        self,
        algorithm=None,
        algorithm_immune=None,
        gen=1,
        select_method=cstrs_immune_system.select_method.BEST_ANTIBODY,
        inject_method=cstrs_immune_system.inject_method.CHAMPION,
        distance_method=cstrs_immune_system.distance_method.EUCLIDEAN,
        phi=0.5,
        gamma=0.5,
        sigma=1. / 3.,
        f_tol=1e-15,
        x_tol=1e-15):
    """
    Constructs an immune system algorithm for constrained optimization.

    USAGE: algorithm._cstrs_immune_system(algorithm = _algorithm.jde(), algorithm_immune = _algorithm.jde(), gen = 1, select_method = cstrs_immune_system.select_method.BEST_ANTIBODY, inject_method = cstrs_immune_system.inject_method.CHAMPION, distance_method = cstrs_immune_system.distance_method.EUCLIDEAN, phi = 0.5, gamma = 0.5, sigma = 1./3., ftol = 1e-15, xtol = 1e-15):

    * algorithm: optimizer to use as 'original' optimization method. Its number of generations should be set to 1.
    * algorithm_2: optimizer to use as 'original' optimization method for the evolution of the immune system.
    * gen: number of generations.
    * select_method: cstrs_immune_system.select_method.BEST_ANTIBODY by default, the method used for selecting the antibodies.
    * inject_method: cstrs_immune_system.inject_method.CHAMPION by default, the method used for reinjecting the antibodies.
    * distance_method: cstrs_immune_system.distance_method.EUCLIDEAN by default, the method used for computing the distance to the antigenes population.
    * Two possibilities are available: CHAMPION, and BEST25.
    * phi: 0.5 by default. The feasible fraction selection to compute the mean value.
    * gamma: 0.5 by default. The number of antigens selected / number of total antigens.
    * sigma: 1/3 by default. The number of antibodies / number of antigens.
    * ftol: 1e-15 by default. The stopping criteria on the x tolerance.
    * xtol: 1e-15 by default. The stopping criteria on the f tolerance.
    """
    arg_list = []

    if algorithm is None:
        algorithm = _algorithm.jde()
    if algorithm_immune is None:
        algorithm_immune = _algorithm.jde()
    arg_list.append(algorithm)
    arg_list.append(algorithm_immune)
    arg_list.append(gen)
    arg_list.append(select_method)
    arg_list.append(inject_method)
    arg_list.append(distance_method)
    arg_list.append(phi)
    arg_list.append(gamma)
    arg_list.append(sigma)
    arg_list.append(f_tol)
    arg_list.append(x_tol)
    self._orig_init(*arg_list)
cstrs_immune_system._orig_init = cstrs_immune_system.__init__
cstrs_immune_system.__init__ = _cstrs_immune_system_ctor


def _cstrs_core_ctor(
        self,
        algorithm=None,
        repair_algorithm=None,
        gen=1,
        repair_frequency=10,
        repair_ratio=1.,
        f_tol=1e-15,
        x_tol=1e-15):
    """
    Constructs CORE (Constrained Optimization by Random Evolution) algorithm for constrained optimization (belong to the family of repairing techniques).

    USAGE: algorithm._cstrs_core(algorithm = _algorithm.jde(), repair_algorithm = _algorithm.jde(), gen = 1, repair_frequency = 10, repair_ratio = 1., f_tol = 1e-15, x_tol = 1e-15):

    * algorithm: optimizer to use as 'original' optimization method. Its number of generations should be set to 1.
    * repair_algorithm: optimizer to use as 'repairing' algorithm. It should be able to deal with population of size 1.
    * gen: number of generations.
    * repair_frequency: The infeasible are repaired at each repair frequency generations.
    * repair_ratio: ratio of repaired individuals over infeasible (a ratio of 1 will repair all the individuals).
    * ftol: 1e-15 by default. The stopping criteria on the x tolerance.
    * xtol: 1e-15 by default. The stopping criteria on the f tolerance.
    """
    arg_list = []
    if algorithm is None:
        algorithm = _algorithm.jde()
    if repair_algorithm is None:
        repair_algorithm = _algorithm.jde()
    arg_list.append(algorithm)
    arg_list.append(repair_algorithm)
    arg_list.append(gen)
    arg_list.append(repair_frequency)
    arg_list.append(repair_ratio)
    arg_list.append(f_tol)
    arg_list.append(x_tol)
    self._orig_init(*arg_list)
cstrs_core._orig_init = cstrs_core.__init__
cstrs_core.__init__ = _cstrs_core_ctor


def _ihs_ctor(
        self,
        iter=100,
        hmcr=0.85,
        par_min=0.35,
        par_max=0.99,
        bw_min=1E-5,
        bw_max=1):
    """
    Constructs an Improved Harmony Search Algorithm

    USAGE: algorithm.ihs(iter = 100, hmcr = 0.85, par_min = 0.35, par_max = 0.99, bw_min = 1E-5, bw_max = 1);

    * iter: number of iterations (improvisations)
    * hmcr: rate of choosing from memory (in ]0,1[)
    * par_min: minimum pitch adjustment rate (in ]0,1[)
    * par_max: maximum pitch adjustment rate (in ]0,1[, > par_min)
    * bw_min: minimum distance bandwidth
    * bw_max: maximum distance bandwidth (> bw_min)
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(iter)
    arg_list.append(hmcr)
    arg_list.append(par_min)
    arg_list.append(par_max)
    arg_list.append(bw_min)
    arg_list.append(bw_max)
    self._orig_init(*arg_list)
ihs._orig_init = ihs.__init__
ihs.__init__ = _ihs_ctor


def _cmaes_ctor(
        self,
        gen=500,
        cc=-1,
        cs=-1,
        c1=-1,
        cmu=-1,
        sigma0=0.5,
        ftol=1e-6,
        xtol=1e-6,
        memory=False,
        screen_output=False):
    """
    Constructs a Covariance Matrix Adaptation Evolutionary Strategy (C++)

    USAGE: algorithm.cmaes(gen = 500, cc = -1, cs = -1, c1 = -1, cmu = -1, sigma0=0.5, ftol = 1e-6, xtol = 1e-6, memory = False, screen_output = False)

    NOTE: In our variant of the algorithm, particle memory is used to extract the elite and reinsertion
    is made aggressively ..... getting rid of the worst guy). Also, the bounds of the problem
    are enforced, as to allow PaGMO machinery to work. Fine control on each iteration can be achieved
    by calling the algo with memory=True and gen=1

    * gen: number of generations
    * cc: time constant for C cumulation (in [0,1]) if -1 automatic values are set
    * cs: time constant for sigma cumulation (in [0,1]) if -1 automatic values are set
    * c1: learning rate for rank-1 update (in [0,1]) if -1 automatic values are set
    * cmu: learning rate for rank-mu update (in [0,1]) if -1 automatic values are set
    * sigma0: starting step (std)
    * xtol: stopping criteria on the x tolerance
    * ftol: stopping criteria on the f tolerance
    * memory: if True the algorithm internal state is saved and used for the next call
    * screen_output: activates screen output of the algorithm (do not use in archipealgo, otherwise the screen will be flooded with
    * 		 different island outputs)
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(gen)
    arg_list.append(cc)
    arg_list.append(cs)
    arg_list.append(c1)
    arg_list.append(cmu)
    arg_list.append(sigma0)
    arg_list.append(ftol)
    arg_list.append(xtol)
    arg_list.append(memory)
    self._orig_init(*arg_list)
    self.screen_output = screen_output
cmaes._orig_init = cmaes.__init__
cmaes.__init__ = _cmaes_ctor


def _inverover_ctor(self, gen=100000, ri=0.05, type="random"):
    """
    Constructs a Inverover algorithm

    REF: "Inver-over Operator for the TSP"
    G Tao, Z Michalewicz, Parallel Problem Solving from Nature - PPSN V, 1998.
    https://cs.adelaide.edu.au/~zbyszek/Papers/p44.pdf

    USAGE: algorithm.inverover(gen=100000, ri=0.05, type="random")

    * gen: number of generations
    * ri: probability for a random inversion (mutation probability)
    * ini_type: algorithm that is used for the initialization of the population
           1. "random"	random initialization with feasible tours
           2. "nn"	using the Nearest-Neighbor algorithm
    """

    from PyGMO.algorithm._algorithm import _tsp_ini_type

    def initialization_type(x):
        return {
            "random": _tsp_ini_type.RANDOM,
            "nn": _tsp_ini_type.NN
        }[x]

    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(gen)
    arg_list.append(ri)
    arg_list.append(initialization_type(type))
    self._orig_init(*arg_list)
inverover._orig_init = inverover.__init__
inverover.__init__ = _inverover_ctor


def _monte_carlo_ctor(self, iter=10000):
    """
    Constructs a Monte Carlo Algorithm

    USAGE: algorithm.monte_carlo(iter = 10000)

    NOTE: At the end of each iteration, the randomly generated
            point substitutes the worst in the population if better

    * iter: number of Monte Carlo runs
    """
    # We set the defaults or the kwargs
    arg_list = []
    arg_list.append(iter)
    self._orig_init(*arg_list)
monte_carlo._orig_init = monte_carlo.__init__
monte_carlo.__init__ = _monte_carlo_ctor

# NLOPT algorithms (only if PyGMO has been compiled woth nlopt option
# activated)
if "nlopt" in str(_get_algorithm_list()):
    def _nlopt_bobyqa_ctor(self, max_iter=100, ftol=1e-6, xtol=1e-6):
        """
        Constructs a BOBYQA algorithm (Bound Optimization BY Quadratic Approximation) (NLOPT)

        USAGE: algorithm.nlopt_bobyqa(max_iter = 100, ftol = 1e-6, xtol = 1e-6)

        * max_iter: stop-criteria (number of iterations)
        * ftol: stop-criteria (absolute on the obj-fun)
        * xtol: stop-criteria (absolute on the chromosome)
        """
        # We set the defaults or the kwargs
        arg_list = []
        arg_list.append(max_iter)
        arg_list.append(ftol)
        arg_list.append(xtol)
        self._orig_init(*arg_list)
    nlopt_bobyqa._orig_init = nlopt_bobyqa.__init__
    nlopt_bobyqa.__init__ = _nlopt_bobyqa_ctor

    def _nlopt_sbplx_ctor(self, max_iter=100, ftol=1e-6, xtol=1e-6):
        """
        Constructs a Subplex (a variant of Nelder-Mead that uses Nelder-Mead on a sequence of subspaces) (NLOPT)

        USAGE: algorithm.nlopt_sbplx(max_iter = 100, ftol = 1e-6, xtol = 1e-6)

        * max_iter: stop-criteria (number of iterations)
        * ftol: stop-criteria (absolute on the obj-fun)
        * xtol: stop-criteria (absolute on the chromosome)
        """
        # We set the defaults or the kwargs
        arg_list = []
        arg_list.append(max_iter)
        arg_list.append(ftol)
        arg_list.append(xtol)
        self._orig_init(*arg_list)
    nlopt_sbplx._orig_init = nlopt_sbplx.__init__
    nlopt_sbplx.__init__ = _nlopt_sbplx_ctor

    def _nlopt_cobyla_ctor(self, max_iter=100, ftol=1e-6, xtol=1e-6):
        """
        Constructs a Constrained Optimization BY Linear Approximation (COBYLA) algorithm (NLOPT)

        USAGE: algorithm.nlopt_cobyla(max_iter = 100, ftol = 1e-6, xtol = 1e-6)

        * max_iter: stop-criteria (number of iterations)
        * ftol: stop-criteria (absolute on the obj-fun)
        * xtol: stop-criteria (absolute on the chromosome)
        """
        # We set the defaults or the kwargs
        arg_list = []
        arg_list.append(max_iter)
        arg_list.append(ftol)
        arg_list.append(xtol)
        self._orig_init(*arg_list)
    nlopt_cobyla._orig_init = nlopt_cobyla.__init__
    nlopt_cobyla.__init__ = _nlopt_cobyla_ctor

    def _nlopt_mma_ctor(self, max_iter=100, ftol=1e-6, xtol=1e-6):
        """
        Constructs a Method of Moving Asymptotes (MMA) algorithm (NLOPT)

        USAGE: algorithm.nlopt_mma(max_iter = 100, ftol = 1e-6, xtol = 1e-6)

        * max_iter: stop-criteria (number of iterations)
        * ftol: stop-criteria (absolute on the obj-fun)
        * xtol: stop-criteria (absolute on the chromosome)
        """
        # We set the defaults or the kwargs
        arg_list = []
        arg_list.append(max_iter)
        arg_list.append(ftol)
        arg_list.append(xtol)
        self._orig_init(*arg_list)
    nlopt_mma._orig_init = nlopt_mma.__init__
    nlopt_mma.__init__ = _nlopt_mma_ctor

    def _nlopt_auglag_ctor(
            self,
            aux_algo_id=1,
            max_iter=100,
            ftol=1e-6,
            xtol=1e-6,
            aux_max_iter=100,
            aux_ftol=1e-6,
            aux_xtol=1e-6):
        """
        Constructs an Augmented agrangian Algotihm (NLOPT)

        USAGE: algorithm.nlopt_mma(aux_algo_id = 1, max_iter = 100, ftol = 1e-6, xtol = 1e-6, aux_max_iter = 100, aux_ftol = 1e-6, aux_xtol = 1e-6)

        * aux_algo_id: auxiliary  optimizer id
                1: SBPLX
                2: COBYLA
                3: BOBYQA
                4: Low Storage BFGS
        * max_iter: stop-criteria (number of iterations)
        * ftol: stop-criteria (absolute on the obj-fun)
        * xtol: stop-criteria (absolute on the chromosome)
        * aux_max_iter: stop-criteria for the auxiliary optimizer (number of iterations)
        * aux_ftol: stop-criteria for the auxiliary optimizer (absolute on the obj-fun)
        * aux_xtol: stop-criteria for the auxiliary optimizer (absolute on the chromosome)
        """
        # We set the defaults or the kwargs
        arg_list = []
        arg_list.append(aux_algo_id)
        arg_list.append(max_iter)
        arg_list.append(ftol)
        arg_list.append(xtol)
        arg_list.append(aux_max_iter)
        arg_list.append(aux_ftol)
        arg_list.append(aux_xtol)
        self._orig_init(*arg_list)
    nlopt_auglag._orig_init = nlopt_auglag.__init__
    nlopt_auglag.__init__ = _nlopt_auglag_ctor

    def _nlopt_auglag_eq_ctor(
            self,
            aux_algo_id=1,
            max_iter=100,
            ftol=1e-6,
            xtol=1e-6,
            aux_max_iter=100,
            aux_ftol=1e-6,
            aux_xtol=1e-6):
        """
        Constructs an Augmented agrangian Algotihm (using penalties only for the equalities) (NLOPT)

        USAGE: algorithm.nlopt_auglag_eq(aux_algo_id = 1, max_iter = 100, ftol = 1e-6, xtol = 1e-6, aux_max_iter = 100, aux_ftol = 1e-6, aux_xtol = 1e-6)

        * aux_algo_id: auxiliary (local) optimizer id
                1: COBYLA
                2: MMA
        * max_iter: stop-criteria (number of iterations)
        * ftol: stop-criteria (absolute on the obj-fun)
        * xtol: stop-criteria (absolute on the chromosome)
        * aux_max_iter: stop-criteria for the auxiliary optimizer (number of iterations)
        * aux_ftol: stop-criteria for the auxiliary optimizer (absolute on the obj-fun)
        * aux_xtol: stop-criteria for the auxiliary optimizer (absolute on the chromosome)
        """
        # We set the defaults or the kwargs
        arg_list = []
        arg_list.append(aux_algo_id)
        arg_list.append(max_iter)
        arg_list.append(ftol)
        arg_list.append(xtol)
        arg_list.append(aux_max_iter)
        arg_list.append(aux_ftol)
        arg_list.append(aux_xtol)
        self._orig_init(*arg_list)
    nlopt_auglag_eq._orig_init = nlopt_auglag_eq.__init__
    nlopt_auglag_eq.__init__ = _nlopt_auglag_eq_ctor

    def _nlopt_slsqp_ctor(self, max_iter=100, ftol=1e-6, xtol=1e-6):
        """
        Constructs a Sequential Least SQuares Programming algorithm (SLSQP) algorithm (NLOPT)

        USAGE: algorithm.nlopt_slsqp(max_iter = 100, ftol = 1e-6, xtol = 1e-6)

        * max_iter: stop-criteria (number of iterations)
        * ftol: stop-criteria (absolute on the obj-fun)
        * xtol: stop-criteria (absolute on the chromosome)
        """
        # We set the defaults or the kwargs
        arg_list = []
        arg_list.append(max_iter)
        arg_list.append(ftol)
        arg_list.append(xtol)
        self._orig_init(*arg_list)
    nlopt_slsqp._orig_init = nlopt_slsqp.__init__
    nlopt_slsqp.__init__ = _nlopt_slsqp_ctor

# GSL algorithms (only if PyGMO has been compiled with gsl option activated)
if "gsl" in str(_get_algorithm_list()):
    def _gsl_bfgs_ctor(
            self,
            max_iter=100,
            step_size=1e-8,
            tol=1e-8,
            grad_step_size=0.01,
            grad_tol=0.0001):
        """
        Constructs a BFGS Algorithm (GSL)

        USAGE: algorithm.gsl_bfgs(max_iter = 100, step_size = 1e-8, tol = 1e-8, grad_step_size = 0.01, grad_tol = 0.0001)

        * max_iter: maximum number of iterations
        * step_size: size of the first trial step.
        * tol: accuracy of the line minimisation.
        * grad_step_size: step size for the numerical computation of the gradient.
        * grad_tol: tolerance when testing the norm of the gradient as stopping criterion.
        """
        # We set the defaults or the kwargs
        arg_list = []
        arg_list.append(max_iter)
        arg_list.append(grad_tol)
        arg_list.append(grad_step_size)
        arg_list.append(tol)
        arg_list.append(step_size)
        self._orig_init(*arg_list)
    gsl_bfgs._orig_init = gsl_bfgs.__init__
    gsl_bfgs.__init__ = _gsl_bfgs_ctor

    def _gsl_bfgs2_ctor(
            self,
            max_iter=100,
            step_size=1e-8,
            tol=1e-8,
            grad_step_size=0.01,
            grad_tol=0.0001):
        """
        Constructs a BFGS2 Algorithm (GSL)

        NOTE: in GSL, BFGS2 is a more efficient version of BFGS

        USAGE: algorithm.gsl_bfgs2(max_iter = 100, step_size = 1e-8, tol = 1e-8, grad_step_size = 0.01, grad_tol = 0.0001);

        * max_iter: maximum number of iterations
        * step_size: size of the first trial step.
        * tol: accuracy of the line minimisation.
        * grad_step_size: step size for the numerical computation of the gradient.
        * grad_tol: tolerance when testing the norm of the gradient as stopping criterion.
        """
        # We set the defaults or the kwargs
        arg_list = []
        arg_list.append(max_iter)
        arg_list.append(grad_tol)
        arg_list.append(grad_step_size)
        arg_list.append(tol)
        arg_list.append(step_size)
        self._orig_init(*arg_list)
    gsl_bfgs2._orig_init = gsl_bfgs2.__init__
    gsl_bfgs2.__init__ = _gsl_bfgs2_ctor

    def _gsl_fr_ctor(
            self,
            max_iter=100,
            step_size=1e-8,
            tol=1e-8,
            grad_step_size=0.01,
            grad_tol=0.0001):
        """
        Constructs a Fletcher-Reeves conjugate gradient (GSL)

        USAGE: algorithm.gsl_fr(max_iter = 100, step_size = 1e-8, tol = 1e-8, grad_step_size = 0.01, grad_tol = 0.0001)

        * max_iter: maximum number of iterations
        * step_size: size of the first trial step.
        * tol: accuracy of the line minimisation.
        * grad_step_size: step size for the numerical computation of the gradient.
        * grad_tol: tolerance when testing the norm of the gradient as stopping criterion.
        """
        # We set the defaults or the kwargs
        arg_list = []
        arg_list.append(max_iter)
        arg_list.append(grad_tol)
        arg_list.append(grad_step_size)
        arg_list.append(tol)
        arg_list.append(step_size)
        self._orig_init(*arg_list)
    gsl_fr._orig_init = gsl_fr.__init__
    gsl_fr.__init__ = _gsl_fr_ctor

    def _gsl_pr_ctor(
            self,
            max_iter=100,
            step_size=1e-8,
            tol=1e-8,
            grad_step_size=0.01,
            grad_tol=0.0001):
        """
        Constructs a Polak-Ribiere conjugate gradient (GSL)

        USAGE: algorithm.gsl_pr2(max_iter = 100, step_size = 1e-8, tol = 1e-8, grad_step_size = 0.01, grad_tol = 0.0001);

        * max_iter: maximum number of iterations
        * step_size: size of the first trial step.
        * tol: accuracy of the line minimisation.
        * grad_step_size: step size for the numerical computation of the gradient.
        * grad_tol: tolerance when testing the norm of the gradient as stopping criterion.
        """
        # We set the defaults or the kwargs
        arg_list = []
        arg_list.append(max_iter)
        arg_list.append(grad_tol)
        arg_list.append(grad_step_size)
        arg_list.append(tol)
        arg_list.append(step_size)
        self._orig_init(*arg_list)
    gsl_pr._orig_init = gsl_pr.__init__
    gsl_pr.__init__ = _gsl_pr_ctor

    def _gsl_nm_ctor(self, max_iter=100, step_size=1e-8, tol=1e-8):
        """
        Constructs a Nelder-Mead Algorithm (GSL)

        USAGE: algorithm.gsl_nm(max_iter = 100, step_size = 1e-8, tol = 1e-8);

        * max_iter: maximum number of iterations
        * step_size: size of the first trial step.
        * tol: accuracy of the line minimisation.
        """
        # We set the defaults or the kwargs
        arg_list = []
        arg_list.append(max_iter)
        arg_list.append(tol)
        arg_list.append(step_size)
        self._orig_init(*arg_list)
    gsl_nm._orig_init = gsl_nm.__init__
    gsl_nm.__init__ = _gsl_nm_ctor

    def _gsl_nm2_ctor(self, max_iter=100, step_size=1e-8, tol=1e-8):
        """
        Constructs a Nelder-Mead algorithm (Variant2) (GSL)

        USAGE: algorithm.gsl_nm2(max_iter = 100, step_size = 1e-8, tol = 1e-8)

        * max_iter: maximum number of iterations
        * step_size: size of the first trial step.
        * tol: accuracy of the line minimisation.
        """
        # We set the defaults or the kwargs
        arg_list = []
        arg_list.append(max_iter)
        arg_list.append(tol)
        arg_list.append(step_size)
        self._orig_init(*arg_list)
    gsl_nm2._orig_init = gsl_nm2.__init__
    gsl_nm2.__init__ = _gsl_nm2_ctor

    def _gsl_nm2rand_ctor(self, max_iter=100, step_size=1e-8, tol=1e-8):
        """
        Constructs a Nelder-Mead algorithm (Variant2 + randomly oriented initial simplex) (GSL)

        USAGE: algorithm.gsl_nm2rand(max_iter = 100, step_size = 1e-8, tol = 1e-8);

        * max_iter: maximum number of iterations
        * step_size: size of the first trial step.
        * tol: accuracy of the line minimisation.
        """
        # We set the defaults or the kwargs
        arg_list = []
        arg_list.append(max_iter)
        arg_list.append(tol)
        arg_list.append(step_size)
        self._orig_init(*arg_list)
    gsl_nm2rand._orig_init = gsl_nm2rand.__init__
    gsl_nm2rand.__init__ = _gsl_nm2rand_ctor

# IPOPT algorithm (only if PyGMO has been compiled with the ipopt option
# activated)
if "ipopt" in str(_get_algorithm_list()):
    def _ipopt_ctor(
            self,
            max_iter=100,
            constr_viol_tol=1e-08,
            dual_inf_tol=1e-08,
            compl_inf_tol=1e-08,
            nlp_scaling_method=True,
            obj_scaling_factor=1.0,
            mu_init=0.1,
            screen_output=False):
        """
        Constructs an Interior Point OPTimization Algorithm (IPOPT)

        USAGE: algorithm.ipopt(major_iter = 100, constr_viol_tol = 1e-08, dual_inf_tol = 1e-08, compl_inf_tol = 1e-08, screen_output = False);

        * max_iter: Maximum number of major iterations
        * constr_viol_tol: Constraint violation tolerance
        * dual_inf_tol: Dual infeasibility tolerance
        * compl_inf_tol: Complementary feasibility tolerance
        * nlp_scaling_method Select if the "gradient-based" scaling of the  NLP should be used
        * obj_scaling_factor Scaling factor for the objective function.
        * mu_init Initial value for the barrier parameter.
        * screen_output: Activates output on screen
        """
        # We set the defaults or the kwargs
        arg_list = []
        arg_list.append(max_iter)
        arg_list.append(constr_viol_tol)
        arg_list.append(dual_inf_tol)
        arg_list.append(compl_inf_tol)
        arg_list.append(nlp_scaling_method)
        arg_list.append(obj_scaling_factor)
        arg_list.append(mu_init)
        self._orig_init(*arg_list)
        self.screen_output = screen_output
    ipopt._orig_init = ipopt.__init__
    ipopt.__init__ = _ipopt_ctor

# SNOPT algorithm (only if PyGMO has been compiled with the snopt option
# activated)
if "snopt" in str(_get_algorithm_list()):
    def _snopt_ctor(
            self,
            major_iter=100,
            feas_tol=1e-8,
            opt_tol=1e-6,
            screen_output=False):
        """
        Constructs SNOPT Algorithm

        USAGE: algorithm.snopt(major_iter = 100, feas_tol = 1e-6, opt_tol = 1e-6, screen_output = False);

        * major_iter: Maximum number of major iterations
        * feas_tol: Feasibility tolerance
        * opt_tol: Optimality tolerance
        * screen_output: Activates output on screen
        """
        # We set the defaults or the kwargs
        arg_list = []
        arg_list.append(major_iter)
        arg_list.append(feas_tol)
        arg_list.append(opt_tol)
        self._orig_init(*arg_list)
        self.screen_output = screen_output
    snopt._orig_init = snopt.__init__
    snopt.__init__ = _snopt_ctor

# WORHP algorithm (only if PyGMO has been compiled with the worhp option
# activated)
if "worhp" in str(_get_algorithm_list()):
    def _worhp_ctor(
            self,
            MaxIter=100,
            TolFeas=1e-8,
            TolOpti=1e-6,
            screen_output=False):
        """
        Constructs WORHP Algorithm

        USAGE: algorithm.worhp(screen_output = False);

        * MaxIter: Maximum number of major iterations
        * TolFeas: Feasibility tolerance
        * TolOpti: Optimality tolerance
        * screen_output: Activates output on screen
        """
        # We set the defaults or the kwargs
        arg_list = []
        arg_list.append(MaxIter)
        arg_list.append(TolFeas)
        arg_list.append(TolOpti)
        arg_list.append(screen_output)
        self._orig_init(*arg_list)
    worhp._orig_init = worhp.__init__
    worhp.__init__ = _worhp_ctor
