from PyGMO.problem._problem import zdt, dtlz


def _mo3d_plot(self, pop, a=40, comp=[0, 1, 2]):
    """
    Generic plot-method for multi-objective optimization problems with more then 2 objectives

    USAGE: prob.plot(pop, comp[0,2,3])
    * pop: population of solutions to the problem
    * a: angle of view on which the 3d-plot is created
    * comp: indexes the fitness dimension for x,y and z axis in that order
    """
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt
    import numpy as np

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    fit = np.transpose([ind.cur_f for ind in pop])
    try:
        ax.plot(fit[comp[0]], fit[comp[1]], fit[comp[2]], 'ro')
    except IndexError:
        print('Error. Please choose correct fitness dimensions for printing!')

    ax.view_init(azim=a)
    return ax


def _dtlz234_plot(pop, a=40, comp=[0, 1, 2]):
    """
    Specific plot-method for the DTLZ2, DTLZ3 and DTLZ4 - plotting also the optimal pareto-front

    USAGE: prob.plot(pop, comp[0,2,3])

    * pop: population of solutions to the problem

    * a: angle of view on which the 3d-plot is created

    * comp: indexes the fitness dimension for x,y and z axis in that order
    """

    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt
    import numpy as np

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # plot the wireframe of the known optimal pareto front
    thetas = np.linspace(0, (np.pi / 2.0), 30)
    # gammas = np.linspace(-np.pi / 4, np.pi / 4, 30)
    gammas = np.linspace(0, (np.pi / 2.0), 30)

    x_frame = np.outer(np.cos(thetas), np.cos(gammas))
    y_frame = np.outer(np.cos(thetas), np.sin(gammas))
    z_frame = np.outer(np.sin(thetas), np.ones(np.size(gammas)))

    ax.view_init(azim=a)

    ax.set_autoscalex_on(False)
    ax.set_autoscaley_on(False)
    ax.set_autoscalez_on(False)

    ax.set_xlim(0, 1.8)
    ax.set_ylim(0, 1.8)
    ax.set_zlim(0, 1.8)

    ax.plot_wireframe(x_frame, y_frame, z_frame)

    # plot the individuals of the population
    fit = np.transpose([ind.cur_f for ind in pop])
    try:
        ax.plot(fit[comp[0]], fit[comp[1]], fit[comp[2]], 'ro')
    except IndexError:
        print('Error. Please choose correct fitness dimensions for printing!')
    return ax


def _zdt_ctor(self, prob_id=1, param_1=None):
    """
    Constructs a multi-objective box-constrained problem from the ZDT testsuite

    NOTE: K Deb, A Pratap, S Agarwal: A fast and elitist multiobjective genetic algorithm: NSGA-II, IEEE Transactions on, 2002

    USAGE: problem.zdt(prob_id = 1, param_1 = 30)

    * prob_id: Problem number, one of [1,2,...6]
    * param_1: problem dimension for all ZDT problems except ZDT5 (here it is the number of binary strings used)
    """

    arg_list = []
    arg_list.append(prob_id)
    if param_1 is None:
        if prob_id in [1, 2, 3, 4]:
            arg_list.append(30)
        elif prob_id == 5:
            arg_list.append(11)
        else:
            arg_list.append(10)
    else:
        arg_list.append(param_1)
    self._orig_init(*arg_list)

zdt._orig_init = zdt.__init__
zdt.__init__ = _zdt_ctor


def _dtlz_ctor(self, prob_id=1, k=None, fdim=3, alpha=100):
    """
    Constructs a multi-objective box-constrained problem from the DTLZ testsuite

NOTE: K Deb, L Thiele, M Laumanns, E Zitzler, Scalable test problems for evolutionary multiobjective optimization

    USAGE: problem.dtlz(prob_id = 1, k = 20, fdim = 4)

    * prob_id: Problem number, one of [1,2,...7]
    * k: paramter defining integer dimension of the problem: k + fdim - 1
    * fdim: number of objectives
    * alpha: controls density of solutions (just used for prob_id = 4)
    """

    arg_list = []
    arg_list.append(prob_id)
    if k is None:
        if prob_id == 1:
            arg_list.append(5)
        elif prob_id in [2, 3, 4, 5, 6]:
            arg_list.append(10)
        else:
            arg_list.append(20)
    else:
        arg_list.append(k)
    arg_list.append(fdim)
    arg_list.append(alpha)
    self._orig_init(*arg_list)
    if prob_id in [2, 3, 4]:
        self.plot = _dtlz234_plot

dtlz.plot = _mo3d_plot
dtlz._orig_init = dtlz.__init__
dtlz.__init__ = _dtlz_ctor
