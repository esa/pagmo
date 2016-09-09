# -*- coding: utf-8 -*-
from PyGMO.core._core import *
import threading as _threading
import signal as _signal
import os as _os

__doc__ = 'PyGMO core module.'
__all__ = [
    'archipelago',
    'base_island',
    'champion',
    'distribution_type',
    'individual',
    'ipy_island',
    'island',
    'local_island',
    'migration_direction',
    'population',
    'py_island']

_orig_signal = _signal.getsignal(_signal.SIGINT)
_main_pid = _os.getpid()

# Alternative signal handler which ignores sigint if called from a child
# process.


def _sigint_handler(signum, frame):
    import os
    if os.getpid() == _main_pid:
        _orig_signal(signum, frame)

_signal.signal(_signal.SIGINT, _sigint_handler)

# Global lock used when starting processes.
_process_lock = _threading.Lock()

# Raw C++ base island class.
_base_island = _core._base_island


class base_island(_core._base_island):

    def __init__(self, *args):
        if len(args) == 0:
            raise ValueError(
                "Cannot initialise base island without parameters for the constructor.")
        _core._base_island.__init__(self, *args)

    def get_name(self):
        return str(type(self))

    def __get_deepcopy__(self):
        from copy import deepcopy
        return deepcopy(self)


def _generic_island_ctor(self, *args, **kwargs):
    """Unnamed arguments:

            #. algorithm
            #. problem or population
            #. number of individuals (optional and valid only if the second argument is a problem, defaults to 0 if not specified)

    Keyword arguments:

            * *s_policy* -- migration selection policy (defaults to 'best selection' policy)
            * *r_policy* -- migration replacement policy (defaults to 'fair replacement' policy)

    """
    from PyGMO.algorithm._algorithm import _base as _base_algorithm
    from PyGMO.algorithm import base as base_algorithm
    from PyGMO.problem._problem import _base as _base_problem
    from PyGMO.problem._problem import _base_stochastic as _base_problem_stochastic
    from PyGMO.problem import base as base_problem
    from PyGMO.problem import base_stochastic as base_problem_stochastic
    from PyGMO.migration._migration import best_s_policy, fair_r_policy, _base_s_policy, _base_r_policy

    if len(args) < 2 or len(args) > 3:
        raise ValueError(
            "Unnamed arguments list must have either 2 or three elements, but %d elements were found instead." %
            (len(args),))
    if not isinstance(args[0], _base_algorithm):
        raise TypeError("The first unnamed argument must be an algorithm.")
    ctor_args = [args[0]]
    if isinstance(args[1], _base_problem) or isinstance(args[1], _base_problem_stochastic):
        ctor_args.append(args[1])
        if len(args) == 3:
            if not isinstance(args[2], int):
                raise TypeError(
                    "Please provide an integer for the number of individuals in the island.")
            ctor_args.append(args[2])
        else:
            ctor_args.append(0)
    elif isinstance(args[1], population):
        if len(args) == 3:
            raise ValueError(
                "When the second unnamed argument is a population, there cannot be a third unnamed argument.")
        ctor_args.append(args[1])
    else:
        raise TypeError(
            "The second unnamed argument must be either a problem or a population.")

    if 's_policy' in kwargs:
        ctor_args.append(kwargs['s_policy'])
    else:
        ctor_args.append(best_s_policy())
    if not isinstance(ctor_args[-1], _base_s_policy):
        raise TypeError("s_policy must be a migration selection policy.")

    if 'r_policy' in kwargs:
        ctor_args.append(kwargs['r_policy'])
    else:
        ctor_args.append(fair_r_policy())
    if not isinstance(ctor_args[-1], _base_r_policy):
        raise TypeError("r_policy must be a migration replacement policy.")

    if isinstance(self, base_island):
        super(type(self), self).__init__(*ctor_args)
    elif isinstance(self, _base_island):
        self.__original_init__(*ctor_args)
    else:
        assert(self is None)
        n_pythonic_items = 0
        if isinstance(args[0], base_algorithm):
            n_pythonic_items += 1
        if isinstance(args[1], base_problem) or isinstance(args[1], base_problem_stochastic):
            n_pythonic_items += 1
        elif isinstance(args[1], population) and (isinstance(args[1].problem, base_problem) or isinstance(args[1], base_problem_stochastic)):
            n_pythonic_items += 1
        if n_pythonic_items > 0:
            return py_island(*args, **kwargs)
        else:
            return local_island(*args, **kwargs)

local_island.__original_init__ = local_island.__init__
local_island.__init__ = _generic_island_ctor

# This is the function that will be called by the separate process
# spawned from py_island.


def _process_target(q, a, p):
    try:
        tmp = a.evolve(p)
        q.put(tmp)
    except BaseException as e:
        q.put(e)


class py_island(base_island):

    """Python island.

    This island will launch evolutions using the multiprocessing module, available since Python 2.6.
    Each evolution is transparently dispatched to a Python interpreter in a separate process.

    """
    __init__ = _generic_island_ctor

    def _perform_evolution(self, algo, pop):
        try:
            import multiprocessing as mp
            q = mp.Queue()
            # Apparently creating/starting processes is _not_ thread safe:
            # http://bugs.python.org/issue1731717
            # http://stackoverflow.com/questions/1359795/error-while-using-multiprocessing-module-in-a-python-daemon
            # Protect with a global lock.
            with _process_lock:
                process = mp.Process(
                    target=_process_target, args=(q, algo, pop))
                process.start()
            retval = q.get()
            with _process_lock:
                process.join()
            if isinstance(retval, BaseException):
                raise retval
            return retval
        except BaseException as e:
            print('Exception caught during evolution:')
            print(e)
            raise RuntimeError()

    def get_name(self):
        return "Python multiprocessing island"

# This is the function that will be called by the task client
# in ipy_island.


def _maptask_target(a, p):
    try:
        return a.evolve(p)
    except BaseException as e:
        return e


class ipy_island(base_island):

    """Parallel IPython island.

    This island will launch evolutions using IPython's MapTask interface. The evolution will be dispatched
    to IPython engines that, depending on the configuration of IPython/ipcluster, can reside either on the
    local machine or on other remote machines.

    See: http://ipython.scipy.org/doc/stable/html/parallel/index.html

    """
    # NOTE: when using an IPython island, on quitting IPython there might be a warning message
    # reporting an exception being ignored. This seems to be a problem in the foolscap library:
    # http://foolscap.lothar.com/trac/ticket/147
    # Hopefully it will be fixed in the next versions of the library.
    __init__ = _generic_island_ctor

    def _perform_evolution(self, algo, pop):
        try:
            from ipyparallel import Client
            
            # Create client
            rc = Client()
            # Create Load-balanced view
            lbview = rc.load_balanced_view()
            # Run the task
            lbview.block = True
            ar = lbview.apply(_maptask_target, args=(algo, pop))
            # Get retval
            retval = ar.get() 

            if isinstance(retval, BaseException):
                raise retval
            return retval
        except BaseException as e:
            print('Exception caught during evolution:')
            print(e)
            raise RuntimeError()

    def get_name(self):
        return "Parallel IPython island"


def island(*args, **kwargs):
    return _generic_island_ctor(None, *args, **kwargs)

island.__doc__ = '\n'.join(['Island factory function.\n\nThis function will return an instance of an island object\nbuilt according to the following rule: ' +
                            'if the arguments include\neither a pythonic problem or a pythonic algorithm, then an instance\nof :class:`py_island` will be returned; ' +
                            'otherwise, an instance of\n:class:`local_island` will be returned.'] +
                           [s.replace('\t', '') for s in _generic_island_ctor.__doc__.split('\n')])
# The following is necessary for Python 2. s remains in the workspace and will cause the error:
# AttributeError: 'module' object has no attribute 's'
# However in Python 3 s is not anylonger in the workspace and del s will
# cause an exception!
if 's' in globals():
    del s


def _get_island_list():
    from PyGMO import core
    names = [n for n in dir(core) if not n.startswith(
        '_') and not n.startswith('base') and n.endswith('_island')]
    try:
        from ipyparallel import Client
    except ImportError:
        names = [n for n in names if n != 'ipy_island']

    return [core.__dict__[n] for n in names]


def _generic_archi_ctor(self, *args, **kwargs):
    """
    Unnamed arguments (optional):

            #. algorithm
            #. problem
            #. number of islands
            #. number individual in the population

    Keyword arguments:

            * *topology* -- migration topology (defaults to unconnected)
            * *distribution_type* -- distribution_type (defaults to distribution_type.point_to_point)
            * *migration_direction* -- migration_direction (defaults to migration_direction.destination)
    """

    from PyGMO import topology, algorithm, problem
    from difflib import get_close_matches

    if not((len(args) == 4) or (len(args) == 0)):
        raise ValueError(
            "Unnamed arguments list, when present, must be of length 4, but %d elements were found instead" %
            (len(args),))

    # Append everything in the same list of constructor arguments
    ctor_args = []
    for i in args:
        ctor_args.append(i)

    # Pop all known keywords out of kwargs and add a default value if not
    # provided
    # unconnected is default
    ctor_args.append(kwargs.pop('topology', topology.unconnected()))
    # point-to-point is default
    ctor_args.append(
        kwargs.pop('distribution_type', distribution_type.point_to_point))
    # destination is default
    ctor_args.append(
        kwargs.pop('migration_direction', migration_direction.destination))

    # Check for unknown keywords
    kwlist = ['topology', 'distribution_type', 'migration_direction']
    if kwargs:
        s = "The following unknown keyworded argument was passed to the construtor: "
        for kw in kwargs:
            s += kw
            spam = get_close_matches(kw, kwlist)
            if spam:
                s += " (Did you mean %s?), " % spam[0]
            else:
                s += ", "

        raise ValueError(s[:-2])

    # Constructs an empty archipelago with no islands using the C++ constructor
    self.__original_init__(*ctor_args[-3:])

    # We now push back the correct island type if required
    if (len(args)) == 4:
        if not isinstance(args[0], algorithm._base):
            raise TypeError("The first unnamed argument must be an algorithm")
        if not (isinstance(args[1], problem._base) or isinstance(args[1], problem._base_stochastic)):
            raise TypeError("The second unnamed argument must be a problem")
        if not isinstance(args[2], int):
            raise TypeError(
                "The third unnamed argument must be an integer (i.e. number of islands)")
        if not isinstance(args[3], int):
            raise TypeError(
                "The fourth unnamed argument must be an integer (i.e. population size)")
        for n in range(args[2]):
            self.push_back(island(args[0], args[1], args[3]))

archipelago.__original_init__ = archipelago.__init__
archipelago.__init__ = _generic_archi_ctor


def _archipelago_draw(
        self,
        layout='spring',
        n_color='fitness',
        n_size=15,
        n_alpha=0.5,
        e_alpha=0.1,
        e_arrows=False,
        scale_by_degree=False,
        cmap='default'):
    """
    Draw a visualization of the archipelago using networkx.

    USAGE: pos = archipelago.draw(layout = 'spring', color = 'fitness', n_size = 15, scale_by_degree = False, n_alpha = 0.5, e_alpha = 0.1, cmap = 'default', e_arrows=False)

    * layout: Network layout. Can be 'spring' or 'circular' or a list of values pos returned
            by a previous call of the method (so that positions of the islands can be kept fixed.
    * n_color = Defines the color code for the nodes. Can be one of 'fitness', 'links', ... or the standard matplotlib 'blue' .. etc.
    * n_size: The size of nodes. Becomes scaling factor when scale_by_degree=True.
    * n_alpha: Transparency of nodes. Takes value between 0 and 1.
    * e_arrows: Plots arrows on the edges for directed graphs
    * e_elpha: Transparency of edges. Takes value between 0 and 1.
    * scale_by_degree: When True, nodes will be sized proportional to their degree.
    * cmap: color map. one in matplotlib.pyplot.cm
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError('Could not import the networkx module.')
    try:
        import matplotlib.pyplot as pl
    except ImportError:
        raise ImportError('Could not improt the MatPlotLib module.')

    # We set the graph in networkx
    t = self.topology
    G = t.to_networkx()

    # We scale the node sizes
    node_sizes = list(range(nx.number_of_nodes(G)))
    for i in range(nx.number_of_nodes(G)):
        if scale_by_degree:
            node_sizes[i] = nx.degree(G, i) * n_size
        else:
            node_sizes[i] = n_size

    # We compute the layout
    if layout == 'spring':
        pos = nx.spring_layout(G)
    elif layout == "circular":
        pos = nx.circular_layout(G)
    else:
        pos = layout

    # We compute the color_code
    if n_color == 'fitness':
        node_colors = [-isl.population.champion.f[0] for isl in self]
        m = min(node_colors)
        M = max(node_colors)
    elif n_color == 'links':
        m = min(node_colors)
        M = max(node_colors)
        node_colors = [
            t.get_num_adjacent_vertices(i) for i in range(len(self))]
    elif n_color == 'rank':
        vec = [-isl.population.champion.f[0] for isl in self]
        node_colors = sorted(list(range(len(vec))), key=vec.__getitem__)
        M = max(node_colors)
        m = min(node_colors)

    else:
        node_colors = n_color
        m = 0
        M = 0

    if not m == M:
        node_colors = [(node_colors[i] - float(m)) / (M - m)
                       for i in range(len(self))]

    # And we draw the archipelago .....
    ax = pl.figure()
    if cmap == 'default':
        cmap = pl.cm.Reds_r
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=list(
            range(
                len(self))),
        node_color=node_colors,
        cmap=cmap,
        node_size=node_sizes,
        alpha=n_alpha)
    nx.draw_networkx_edges(G, pos, alpha=e_alpha, arrows=e_arrows)
    pl.axis('off')
    pl.show()
    return pos
archipelago.draw = _archipelago_draw


def _pop_ctor(self, prob_or_pop, n_individuals=0, seed=None):
    """
    Constructs a population.

    Popopulation can be constructed in two ways, specified by prob_or_pop.  If
    prob_or_pop is a population (see USAGE 2 below), the other two arguments
    are ignored, and the population is constructed by performing a deep-copy of
    the provided population.

    USAGE 1: pop = population(problem.zdt(), n_individuals=100, seed=1234)

    * prob_or_prob: problem to be associated with the population
    * n_individuals: number of individuals in the population
    * seed: seed used to randomly initialize the individuals

    USAGE 2:
            from PyGMO import *
            pop1 = population(problem.schwefel(50), 10) #population with 10 individuals
            pop2 = population(pop1) # pop2 is a copy of pop1
    """

    arg_list = []
    arg_list.append(prob_or_pop)
    # For construction by copying, ignore the rest of the arguments (could be
    # the default kwargs).
    if not isinstance(prob_or_pop, population):
        arg_list.append(n_individuals)
        if seed is not None:
            arg_list.append(seed)
    return self._original_init(*arg_list)
population._original_init = population.__init__
population.__init__ = _pop_ctor


def _pop_plot_pareto_fronts(
        pop,
        rgb=(
            0,
            0,
            0),
        comp = [
            0,
            1],
        symbol = 'o',
        size = 6,
        fronts=[]):
    """
    Plots the population pareto front in a 2-D graph

    USAGE: pop.plot_pareto_front(comp = [0,1], rgb=(0,1,0))

    * comp: components of the fitness function to plot in the 2-D window
    * rgb: specify the color of the 1st front (use strong colors here)
    * symbol: marker for the individual
    * size: size of the markersymbol
    * fronts: list of fronts to be plotted (use [0] to only show the first)
    """
    from numpy import linspace
    import matplotlib.pyplot as plt

    if len(comp) != 2:
        raise ValueError(
            'Invalid components of the objective function selected for plot')

    p_dim = pop.problem.f_dimension

    if p_dim == 1:
        raise ValueError(
            'Pareto fronts of a 1-dimensional problem cannot be plotted')

    if not all([c in range(0, p_dim) for c in comp]):
        raise ValueError(
            'You need to select valid components of the objective function')

    p_list = pop.compute_pareto_fronts()
    if (len(fronts) > 0):
        n = len(p_list)
        consistent = [d < n for d in fronts]
        if consistent.count(False) > 0:
            raise ValueError(
                'Check your fronts list, there seem to be not enough fronts')
        p_list = [p_list[idx] for idx in fronts]

    cl = list(zip(linspace(0.9 if rgb[0] else 0.1, 0.9, len(p_list)),
                  linspace(0.9 if rgb[1] else 0.1, 0.9, len(p_list)),
                  linspace(0.9 if rgb[2] else 0.1, 0.9, len(p_list))))

    for id_f, f in enumerate(p_list):
        for ind in f:
            plt.plot([pop[ind].cur_f[comp[0]]],
                     [pop[ind].cur_f[comp[1]]],
                     symbol,
                     color=cl[id_f], markersize=size)
        x = [pop[ind].cur_f[comp[0]] for ind in f]
        y = [pop[ind].cur_f[comp[1]] for ind in f]
        tmp = [(a, b) for a, b in zip(x, y)]
        tmp = sorted(tmp, key=lambda k: k[0])
        plt.step([c[0] for c in tmp], [c[1]
                 for c in tmp], color=cl[id_f], where='post')
    return plt.gca()

population.plot_pareto_fronts = _pop_plot_pareto_fronts


def _pop_race(self, n_winners, min_trials=0, max_feval=500,
              delta=0.05, racers_idx=[], race_best=True, screen_output=False):
    """
    Races individuals in a population

    USAGE: pop.race(n_winners, min_trials = 0, max_feval = 500, delta = 0.05, racers_idx = [], race_best=True, screen_output=False)

    * n_winners: number of winners in the race
    * min_trials: minimum amount of evaluations before an individual can stop racing
* max_feval: budget for objective function evaluation
    * delta: Statistical test confidence
    * racers_idx: indices of the individuals in pop to be raced
    * race_best: when True winners are the best, otherwise winners are the worst
    * screen_output: produces some screen output at each iteration of the race
    """
    arg_list = []
    arg_list.append(n_winners)
    arg_list.append(min_trials)
    arg_list.append(max_feval)
    arg_list.append(delta)
    arg_list.append(racers_idx)
    arg_list.append(race_best)
    arg_list.append(screen_output)
    return self._orig_race(*arg_list)

population._orig_race = population.race
population.race = _pop_race


def _pop_repair(self, idx, repair_algorithm):
    """
    Repairs the individual at the given position

    USAGE: pop.repair(idx, repair_algorithm = _algorithm.jde())

    * idx: index of the individual to repair
repair_algorithm: optimizer to use as 'repairing' algorithm. It should be able to deal with population of size 1.
    """
    arg_list = []
    arg_list.append(idx)
    arg_list.append(repair_algorithm)
    return self._orig_repair(*arg_list)

population._orig_repair = population.repair
population.repair = _pop_repair
