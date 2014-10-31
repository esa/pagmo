from PyGMO.problem._problem import cassini_1, gtoc_1, gtoc_2, cassini_2, rosetta, messenger_full, tandem, laplace, sagas, mga_1dsm_alpha, mga_1dsm_tof, mga_incipit, mga_incipit_cstrs, mga_part
from _problem import _gtoc_2_objective

# Redefining the constructors of all problems to obtain good documentation
# and allowing kwargs


def _cassini_1_ctor(self, objectives=1):
    """
    Constructs a Cassini 1 Problem (Box-Constrained Continuous Single-Objective)

    NOTE: This problem (MGA) belongs to the GTOP database [http://www.esa.int/gsp/ACT/inf/op/globopt.htm]
          Its single objective version has a global minimum at 4.9307 [km/s],
          and it is a deceptive problem with a larger minimum at 5.303 [km/s]

    USAGE: problem.cassini_1(objectives = 1)

    * objectives: number of objectives. 1=DV, 2=DV,DT

    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(objectives)
    self._orig_init(*arg_list)
cassini_1._orig_init = cassini_1.__init__
cassini_1.__init__ = _cassini_1_ctor


def _gtoc_1_ctor(self):
    """
    Constructs a GTOC 1 Problem (Box-Constrained Continuous Single-Objective)

    NOTE: This problem (MGA) belongs to the GTOP database [http://www.esa.int/gsp/ACT/inf/op/globopt.htm]

          Best known global minima is at -1,581,950

    USAGE: problem.gtoc_1()

    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    self._orig_init(*arg_list)
gtoc_1._orig_init = gtoc_1.__init__
gtoc_1.__init__ = _gtoc_1_ctor


def _cassini_2_ctor(self):
    """
    Constructs a Cassini 2 Problem (Box-Constrained Continuous Single-Objective)

    NOTE: This problem (MGA-1DSM) belongs to the GTOP database [http://www.esa.int/gsp/ACT/inf/op/globopt.htm]
          It models the same interplanetary trajectory as the cassini_1 problem, but
          in a more accurate fashion, allowing deep space manouvres

          Best known global minimum is at 8.383 [km/s]

    USAGE: problem.cassini_2()

    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    self._orig_init(*arg_list)
cassini_2._orig_init = cassini_2.__init__
cassini_2.__init__ = _cassini_2_ctor


def _rosetta_ctor(self):
    """
    Constructs a Rosetta Problem (Box-Constrained Continuous Single-Objective)

    NOTE: This problem (MGA-1DSM) belongs to the GTOP database [http://www.esa.int/gsp/ACT/inf/op/globopt.htm]

          Best known global minimum is at 1.343 [km/s]

    USAGE: problem.rosetta()

    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    self._orig_init(*arg_list)
rosetta._orig_init = rosetta.__init__
rosetta.__init__ = _rosetta_ctor


def _messenger_full_ctor(self):
    """
    Constructs a Mesenger Full Problem (Box-Constrained Continuous Single-Objective)

    NOTE: This problem (MGA-1DSM) belongs to the GTOP database [http://www.esa.int/gsp/ACT/inf/op/globopt.htm]

          Best known global minimum is at 2.113

    USAGE: problem.messenger_full()

    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    self._orig_init(*arg_list)
messenger_full._orig_init = messenger_full.__init__
messenger_full.__init__ = _messenger_full_ctor


def _tandem_ctor(self, prob_id=7, max_tof=-1):
    """
    Constructs a TandEM Problem (Box-Constrained Continuous Single-Objective)

    NOTE: This problem (MGA-1DSM) belongs to the GTOP database [http://www.esa.int/gsp/ACT/inf/op/globopt.htm]. The objective function is -log(m_final).

    USAGE: problem.tandem(prob_id = 7, max_tof = -1)

    * prob_id: Selects the problem variant (one of 1..25). All problems differ from the fly-by sequence
    * max_tof = Activates a constriants on the maximum time of flight allowed (in years)

    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(prob_id)
    arg_list.append(max_tof)
    self._orig_init(*arg_list)
tandem._orig_init = tandem.__init__
tandem.__init__ = _tandem_ctor


def _laplace_ctor(self, seq=[3, 2, 3, 3, 5]):
    """
    Constructs a EJSM-Laplace Problem (Box-Constrained Continuous Single-Objective)

    NOTE: This problem (MGA-1DSM) is similar to TandEM, but targets Jupiter and the user
          can specify explicitly the planetary fly-by sequence

    USAGE: problem.laplace(seq = [3,2,3,3,5])

    * seq: The planetary sequence. This is a list of ints that represent the planets to visit
           1 - Mercury, 2 - Venus, 3 - Earth, 4 - Mars, 5 - Jupiter, 6 - Saturn. It must start from 3 (Earth)
           and end with 5 (Jupiter)
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(seq)
    self._orig_init(*arg_list)
laplace._orig_init = laplace.__init__
laplace.__init__ = _laplace_ctor


def _sagas_ctor(self):
    """
    Constructs a SAGAS Problem (Box-Constrained Continuous Single-Objective)

    NOTE: This problem (MGA-1DSM) belongs to the GTOP database [http://www.esa.int/gsp/ACT/inf/op/globopt.htm]

    USAGE: problem.sagas()
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    self._orig_init(*arg_list)
sagas._orig_init = sagas.__init__
sagas.__init__ = _sagas_ctor

gtoc_2.obj = _gtoc_2_objective


def _gtoc_2_ctor(self, ast1=815, ast2=300, ast3=110, ast4=47, n_seg=10, objective=gtoc_2.obj.MASS_TIME):
    """
    Constructs a GTOC 2 Problem (Constrained Continuous Single-Objective)

    NOTE: This problem is a quite faithful transcription of the problem used during the GTOC2 competition
        It Transcribe the whole OCP resulting from the low-thrust dynamics into an NLP. As such it is very
          difficult to find feasible solutions. Note that by default the asteroid sequence is the winning one
          from Turin University.

    USAGE: problem.gtoc_2(ast1 = 815, ast2 = 300, ast3 = 110, ast4 = 47, n_seg = 10, objective = gtoc_2.obj.MASS_TIME)

    * ast1 id of the first asteroid to visit (Group 1:   0 - 95)
    * ast2 id of the second asteroid to visit (Group 2:  96 - 271)
    * ast3 id of the third asteroid to visit (Group 3: 272 - 571)
    * ast4 id of the fourth asteroid to visit (Group 4: 572 - 909)
    * n_seg number of segments to be used per leg
    * obj objective function in the enum {MASS,TIME,MASS_TIME}
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(ast1)
    arg_list.append(ast2)
    arg_list.append(ast3)
    arg_list.append(ast4)
    arg_list.append(n_seg)
    arg_list.append(objective)
    self._orig_init(*arg_list)
gtoc_2._orig_init = gtoc_2.__init__
gtoc_2.__init__ = _gtoc_2_ctor


from PyKEP import planet_ss, epoch, planet_js


def _mga_1dsm_alpha_ctor(
        self, seq=[planet_ss('earth'), planet_ss('venus'), planet_ss('earth')],
        t0=[epoch(0), epoch(1000)], tof=[365.25, 5.0 * 365.25], vinf=[0.5,
                                                                      2.5], multi_objective=False, add_vinf_dep=False, add_vinf_arr=True):
    """
    Constructs an mga_1dsm problem (alpha-encoding)

    USAGE: problem.mga_1dsm(seq = [planet_ss('earth'),planet_ss('venus'),planet_ss('earth')], t0 = [epoch(0),epoch(1000)], tof = [365.25,5.0 * 365.25], vinf = [0.5, 2.5], multi_objective = False, add_vinf_dep = False, add_vinf_arr = True)

    * seq: list of PyKEP planets defining the encounter sequence (including the starting planet)
    * t0: list of two epochs defining the launch window
    * tof: list of two floats defining the minimum and maximum allowed mission length (days)
    * vinf: list of two floats defining the minimum and maximum allowed initial hyperbolic velocity at launch (km/sec)
    * multi_objective: when True constructs a multiobjective problem (dv, T)
    * add_vinf_dep: when True the computed Dv includes the initial hyperbolic velocity (at launch)
    * add_vinf_arr: when True the computed Dv includes the final hyperbolic velocity (at arrival)
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(seq)
    arg_list.append(t0[0])
    arg_list.append(t0[1])
    arg_list.append(tof[0])
    arg_list.append(tof[1])
    arg_list.append(vinf[0])
    arg_list.append(vinf[1])
    arg_list.append(multi_objective)
    arg_list.append(add_vinf_dep)
    arg_list.append(add_vinf_arr)
    self._orig_init(*arg_list)
mga_1dsm_alpha._orig_init = mga_1dsm_alpha.__init__
mga_1dsm_alpha.__init__ = _mga_1dsm_alpha_ctor


def _mga_1dsm_tof_ctor(
    self, seq=[
        planet_ss('earth'), planet_ss('venus'), planet_ss('earth')], t0=[
        epoch(0), epoch(1000)], tof=[
        [
            50, 900], [
            50, 900]], vinf=[
        0.5, 2.5], multi_objective=False, add_vinf_dep=False, add_vinf_arr=True):
    """
    Constructs an mga_1dsm problem (tof-encoding)

    USAGE: problem.mga_1dsm(seq = [planet_ss('earth'),planet_ss('venus'),planet_ss('earth')], t0 = [epoch(0),epoch(1000)], tof = [ [50, 900], [50, 900] ], vinf = [0.5, 2.5], multi_objective = False, add_vinf_dep = False, add_vinf_arr = True)

    * seq: list of PyKEP planets defining the encounter sequence (including the starting planet)
    * t0: list of two epochs defining the launch window
    * tof: list of intervals defining the times of flight (days)
    * vinf: list of two floats defining the minimum and maximum allowed initial hyperbolic velocity at launch (km/sec)
    * multi_objective: when True constructs a multiobjective problem (dv, T)
    * add_vinf_dep: when True the computed Dv includes the initial hyperbolic velocity (at launch)
    * add_vinf_arr: when True the computed Dv includes the final hyperbolic velocity (at arrival)
    """

    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(seq)
    arg_list.append(t0[0])
    arg_list.append(t0[1])
    arg_list.append(tof)
    arg_list.append(vinf[0])
    arg_list.append(vinf[1])
    arg_list.append(multi_objective)
    arg_list.append(add_vinf_dep)
    arg_list.append(add_vinf_arr)
    self._orig_init(*arg_list)
mga_1dsm_tof._orig_init = mga_1dsm_tof.__init__
mga_1dsm_tof.__init__ = _mga_1dsm_tof_ctor


def _mga_incipit_ctor(
    self, seq=[
        planet_js('io'), planet_js('io'), planet_js('europa')], t0=[
        epoch(7305.0), epoch(11323.0)], tof=[
        [
            100, 200], [
            3, 200], [
            4, 100]]):
    """
    USAGE: mga_incipit(seq = [planet_js('io'),planet_js('io'),planet_js('europa')], t0 = [epoch(6905.0),epoch(11323.0)], tof = [[100,200],[3,200],[4,100]])

    * seq: list of jupiter moons defining the trajectory incipit
    * t0:  list of two epochs defining the launch window
    * tof: list of n lists containing the lower and upper bounds for the legs flight times (days)
    """
    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(seq)
    arg_list.append(t0[0])
    arg_list.append(t0[1])
    arg_list.append(tof)
    self._orig_init(*arg_list)
mga_incipit._orig_init = mga_incipit.__init__
mga_incipit.__init__ = _mga_incipit_ctor


def _mga_incipit_cstrs_ctor(
    self, seq=[
        planet_js('io'), planet_js('io'), planet_js('europa')], t0=[
        epoch(7305.0), epoch(11323.0)], tof=[
        [
            100, 200], [
            3, 200], [
            4, 100]], Tmax=300.00, Dmin=2.0):
    """
    USAGE: mga_incipit_cstrs(seq = [planet_js('io'),planet_js('io'),planet_js('europa')], t0 = [epoch(6905.0),epoch(11323.0)], tof = [[100,200],[3,200],[4,100]], Tmax = 365.25, Dmin = 0.2)

    * seq: list of jupiter moons defining the trajectory incipit
    * t0:  list of two epochs defining the launch window
    * tof: list of n lists containing the lower and upper bounds for the legs flight times (days)
    """
    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(seq)
    arg_list.append(t0[0])
    arg_list.append(t0[1])
    arg_list.append(tof)
    arg_list.append(Tmax)
    arg_list.append(Dmin)
    self._orig_init(*arg_list)
mga_incipit_cstrs._orig_init = mga_incipit_cstrs.__init__
mga_incipit_cstrs.__init__ = _mga_incipit_cstrs_ctor


def _mga_part_ctor(
    self, seq=[
        planet_js('europa'), planet_js('europa'), planet_js('europa')], tof=[
        [
            5, 50], [
            5, 50]], t0=epoch(11000), v_inf_in=[
        1500.0, 350.0, 145.0]):
    """
    USAGE: mga_part(seq = [planet_js('europa'),planet_js('europa'),planet_js('europa')], tof = [[5,50],[5,50]], t0 = epoch(11000), v_inf_in[1500.0,350.0,145.0])

    * seq: list of jupiter moons defining the trajectory incipit
    * tof: list of n lists containing the lower and upper bounds for the legs flight times (days)
    * t0:  starting epoch
    * v_inf_in: Incoming spacecraft relative velocity
    """
    # We construct the arg list for the original constructor exposed by
    # boost_python
    arg_list = []
    arg_list.append(seq)
    arg_list.append(tof)
    arg_list.append(t0)
    arg_list.append(v_inf_in)
    self._orig_init(*arg_list)
mga_part._orig_init = mga_part.__init__
mga_part.__init__ = _mga_part_ctor

# Plot of the trajectory for an mga_1dsm problem


def _mga_1dsm_alpha_plot(self, x):
    """
    Plots the trajectory represented by the decision vector x
    """
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from PyKEP.orbit_plots import plot_planet, plot_lambert, plot_kepler
    from PyKEP import epoch, propagate_lagrangian, lambert_problem, fb_prop, AU, MU_SUN, DAY2SEC
    from math import pi, acos, cos, sin
    from scipy.linalg import norm

    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(0, 0, 0, color='y')

    seq = self.get_sequence()

    n = (len(seq) - 1)
    # 1 -  we 'decode' the chromosome recording the various times of flight
    # (days) in the list T
    T = list([0] * (n))
    alpha_sum = 0
    for i in range(n):
        T[i] = x[1] * x[6 + 4 * i]
        alpha_sum += x[6 + 4 * i]
    for i in range(n):
        T[i] /= alpha_sum

    # 2 - We compute the epochs and ephemerides of the planetary encounters
    t_P = list([None] * (n + 1))
    r_P = list([None] * (n + 1))
    v_P = list([None] * (n + 1))
    DV = list([None] * (n + 1))

    for i, planet in enumerate(seq):
        t_P[i] = epoch(x[0] + sum(T[0:i]))
        r_P[i], v_P[i] = planet.eph(t_P[i])
        plot_planet(ax, planet, t0=t_P[i], color=(
            0.8, 0.6, 0.8), legend=True, units = AU)

    # 3 - We start with the first leg
    theta = 2 * pi * x[2]
    phi = acos(2 * x[3] - 1) - pi / 2

    Vinfx = x[4] * cos(phi) * cos(theta)
    Vinfy = x[4] * cos(phi) * sin(theta)
    Vinfz = x[4] * sin(phi)

    v0 = [a + b for a, b in zip(v_P[0], [Vinfx, Vinfy, Vinfz])]
    r, v = propagate_lagrangian(
        r_P[0], v0, x[5] * T[0] * DAY2SEC, seq[0].mu_central_body)
    plot_kepler(
        ax,
        r_P[0],
        v0,
        x[5] *
        T[0] *
        DAY2SEC,
        seq[0].mu_central_body,
        N=100,
        color='b',
        legend=False,
        units=AU)

    # Lambert arc to reach seq[1]
    dt = (1 - x[5]) * T[0] * DAY2SEC
    l = lambert_problem(r, r_P[1], dt, seq[0].mu_central_body)
    plot_lambert(ax, l, sol=0, color='r', legend=False, units=AU)
    v_end_l = l.get_v2()[0]
    v_beg_l = l.get_v1()[0]

    # First DSM occuring at time nu1*T1
    DV[0] = norm([a - b for a, b in zip(v_beg_l, v)])

    # 4 - And we proceed with each successive leg
    for i in range(1, n):
        # Fly-by
        v_out = fb_prop(v_end_l,
                        v_P[i],
                        x[8 + (i - 1) * 4] * seq[i].radius,
                        x[7 + (i - 1) * 4],
                        seq[i].mu_self)
        # s/c propagation before the DSM
        r, v = propagate_lagrangian(
            r_P[i], v_out, x[9 + (i - 1) * 4] * T[i] * DAY2SEC, seq[0].
            mu_central_body)
        plot_kepler(ax,
                    r_P[i],
                    v_out,
                    x[9 + (i - 1) * 4] * T[i] * DAY2SEC,
                    seq[0].mu_central_body,
                    N=100,
                    color='b',
                    legend=False,
                    units=AU)
        # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
        dt = (1 - x[9 + (i - 1) * 4]) * T[i] * DAY2SEC
        l = lambert_problem(r, r_P[i + 1], dt, seq[0].mu_central_body)
        plot_lambert(ax, l, sol=0, color='r', legend=False, units=AU)
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]
        # DSM occurring at time nu2*T2
        DV[i] = norm([a - b for a, b in zip(v_beg_l, v)])
    return ax
mga_1dsm_alpha.plot = _mga_1dsm_alpha_plot


# Plot of the trajectory for an mga_1dsm problem
def _mga_1dsm_tof_plot_old(self, x):
    """
    Plots the trajectory represented by the decision vector x
    """
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from PyKEP.orbit_plots import plot_planet, plot_lambert, plot_kepler
    from PyKEP import epoch, propagate_lagrangian, lambert_problem, fb_prop, AU, MU_SUN, DAY2SEC
    from math import pi, acos, cos, sin
    from scipy.linalg import norm

    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(0, 0, 0, color='y')

    seq = self.get_sequence()

    n = (len(seq) - 1)
    # 1 -  we 'decode' the chromosome recording the various times of flight
    # (days) in the list T
    T = x[5::4]

    # 2 - We compute the epochs and ephemerides of the planetary encounters
    t_P = list([None] * (n + 1))
    r_P = list([None] * (n + 1))
    v_P = list([None] * (n + 1))
    DV = list([None] * (n + 1))

    for i, planet in enumerate(seq):
        t_P[i] = epoch(x[0] + sum(T[0:i]))
        r_P[i], v_P[i] = planet.eph(t_P[i])
        plot_planet(ax, planet, t0=t_P[i], color=(
            0.8, 0.6, 0.8), legend=True, units = AU)

    # 3 - We start with the first leg
    theta = 2 * pi * x[1]
    phi = acos(2 * x[2] - 1) - pi / 2

    Vinfx = x[3] * cos(phi) * cos(theta)
    Vinfy = x[3] * cos(phi) * sin(theta)
    Vinfz = x[3] * sin(phi)

    v0 = [a + b for a, b in zip(v_P[0], [Vinfx, Vinfy, Vinfz])]
    r, v = propagate_lagrangian(
        r_P[0], v0, x[4] * T[0] * DAY2SEC, seq[0].mu_central_body)
    plot_kepler(
        ax,
        r_P[0],
        v0,
        x[4] *
        T[0] *
        DAY2SEC,
        seq[0].mu_central_body,
        N=100,
        color='b',
        legend=False,
        units=AU)

    # Lambert arc to reach seq[1]
    dt = (1 - x[4]) * T[0] * DAY2SEC
    l = lambert_problem(r, r_P[1], dt, seq[0].mu_central_body)
    plot_lambert(ax, l, sol=0, color='r', legend=False, units=AU)
    v_end_l = l.get_v2()[0]
    v_beg_l = l.get_v1()[0]

    # First DSM occuring at time nu1*T1
    DV[0] = norm([a - b for a, b in zip(v_beg_l, v)])

    # 4 - And we proceed with each successive leg
    for i in range(1, n):
        # Fly-by
        v_out = fb_prop(v_end_l,
                        v_P[i],
                        x[7 + (i - 1) * 4] * seq[i].radius,
                        x[6 + (i - 1) * 4],
                        seq[i].mu_self)
        # s/c propagation before the DSM
        r, v = propagate_lagrangian(
            r_P[i], v_out, x[8 + (i - 1) * 4] * T[i] * DAY2SEC, seq[0].
            mu_central_body)
        plot_kepler(ax,
                    r_P[i],
                    v_out,
                    x[8 + (i - 1) * 4] * T[i] * DAY2SEC,
                    seq[0].mu_central_body,
                    N=100,
                    color='b',
                    legend=False,
                    units=AU)
        # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
        dt = (1 - x[8 + (i - 1) * 4]) * T[i] * DAY2SEC
        l = lambert_problem(r, r_P[i + 1], dt, seq[0].mu_central_body)
        plot_lambert(ax, l, sol=0, color='r', legend=False, units=AU)
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]
        # DSM occurring at time nu2*T2
        DV[i] = norm([a - b for a, b in zip(v_beg_l, v)])
    return ax
mga_1dsm_tof.plot_old = _mga_1dsm_tof_plot_old

# Plot of the trajectory of an mga_incipit problem


def _mga_incipit_plot_old(self, x, plot_leg_0=False):
    """
    Plots the trajectory represented by the decision vector x

    Example::

      prob.plot(x)
    """
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from PyKEP.orbit_plots import plot_planet, plot_lambert, plot_kepler
    from PyKEP import epoch, propagate_lagrangian, lambert_problem, fb_prop, AU, MU_SUN, DAY2SEC
    from math import pi, acos, cos, sin
    from scipy.linalg import norm

    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d', aspect='equal')
    ax.scatter(0, 0, 0, color='y')

    JR = 71492000.0
    legs = len(x) / 4
    seq = self.get_sequence()
    common_mu = seq[0].mu_central_body

    # 1 -  we 'decode' the chromosome recording the various times of flight
    # (days) in the list T
    T = x[3::4]

    # 2 - We compute the epochs and ephemerides of the planetary encounters
    t_P = list([None] * legs)
    r_P = list([None] * legs)
    v_P = list([None] * legs)
    DV = list([None] * legs)

    for i, planet in enumerate(seq):
        t_P[i] = epoch(x[0] + sum(T[:i + 1]))
        r_P[i], v_P[i] = planet.eph(t_P[i])
        plot_planet(ax, planet, t0=t_P[i], color=(
            0.8, 0.6, 0.8), legend=True, units = JR)

    # 3 - We start with the first leg: a lambert arc
    theta = 2 * pi * x[1]
    phi = acos(2 * x[2] - 1) - pi / 2
    # phi close to zero is in the moon orbit plane injection
    r = [cos(phi) * sin(theta), cos(phi) * cos(theta), sin(phi)]
    r = [JR * 1000 * d for d in r]

    l = lambert_problem(r, r_P[0], T[0] * DAY2SEC, common_mu, False, False)
    if (plot_leg_0):
        plot_lambert(ax, l, sol=0, color='k', legend=False, units=JR, N=500)

    # Lambert arc to reach seq[1]
    v_end_l = l.get_v2()[0]
    v_beg_l = l.get_v1()[0]

    # 4 - And we proceed with each successive leg
    for i in range(1, legs):
        # Fly-by
        v_out = fb_prop(v_end_l,
                        v_P[i - 1],
                        x[1 + 4 * i] * seq[i - 1].radius,
                        x[4 * i],
                        seq[i - 1].mu_self)
        # s/c propagation before the DSM
        r, v = propagate_lagrangian(
            r_P[i - 1], v_out, x[4 * i + 2] * T[i] * DAY2SEC, common_mu)
        plot_kepler(ax,
                    r_P[i - 1],
                    v_out,
                    x[4 * i + 2] * T[i] * DAY2SEC,
                    common_mu,
                    N=500,
                    color='b',
                    legend=False,
                    units=JR)
        # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
        dt = (1 - x[4 * i + 2]) * T[i] * DAY2SEC
        l = lambert_problem(r, r_P[i], dt, common_mu, False, False)
        plot_lambert(ax, l, sol=0, color='r', legend=False, units=JR, N=500)
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]
    plt.show()
    return ax
mga_incipit.plot_old = _mga_incipit_plot_old

# Plot of the trajectory of an mga_part problem


def _mga_part_plot_old(self, x):
    """
    Plots the trajectory represented by the decision vector x

    Example::

      prob.plot(x)
    """
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from PyKEP.orbit_plots import plot_planet, plot_lambert, plot_kepler
    from PyKEP import epoch, propagate_lagrangian, lambert_problem, fb_prop, AU, MU_SUN, DAY2SEC
    from math import pi, acos, cos, sin
    from scipy.linalg import norm

    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d', aspect='equal')
    ax.scatter(0, 0, 0, color='y')

    JR = 71492000.0
    legs = len(x) / 4
    seq = self.get_sequence()
    common_mu = seq[0].mu_central_body
    start_mjd2000 = self.t0.mjd2000

    # 1 -  we 'decode' the chromosome recording the various times of flight
    # (days) in the list T
    T = x[3::4]

    # 2 - We compute the epochs and ephemerides of the planetary encounters
    t_P = list([None] * (legs + 1))
    r_P = list([None] * (legs + 1))
    v_P = list([None] * (legs + 1))

    for i, planet in enumerate(seq):
        t_P[i] = epoch(start_mjd2000 + sum(T[:i]))
        r_P[i], v_P[i] = planet.eph(t_P[i])
        plot_planet(ax, planet, t0=t_P[i], color=(
            0.8, 0.6, 0.8), legend=True, units = JR)

    v_end_l = [a + b for a, b in zip(v_P[0], self.vinf_in)]
    # 4 - And we iterate on the legs
    for i in range(0, legs):
        # Fly-by
        v_out = fb_prop(v_end_l,
                        v_P[i],
                        x[1 + 4 * i] * seq[i - 1].radius,
                        x[4 * i],
                        seq[i].mu_self)
        # s/c propagation before the DSM
        r, v = propagate_lagrangian(
            r_P[i], v_out, x[4 * i + 2] * T[i] * DAY2SEC, common_mu)
        plot_kepler(ax, r_P[i], v_out, x[4 * i + 2] * T[i] * DAY2SEC,
                    common_mu, N=500, color='b', legend=False, units=JR)
        # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
        dt = (1 - x[4 * i + 2]) * T[i] * DAY2SEC
        l = lambert_problem(r, r_P[i + 1], dt, common_mu, False, False)
        plot_lambert(ax, l, sol=0, color='r', legend=False, units=JR, N=500)
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]
    plt.show()
    return ax
mga_part.plot_old = _mga_part_plot_old

# Plot of concatenated fly-by legs


def _part_plot(x, units, axis, seq, start_mjd2000, vinf_in):
    """
    Plots the trajectory represented by a decision vector x = [beta,rp,eta,T] * N
    associated to a sequence seq, a start_mjd2000 and an incoming vinf_in
    """
    from PyKEP.orbit_plots import plot_planet, plot_lambert, plot_kepler
    from PyKEP import epoch, propagate_lagrangian, lambert_problem, fb_prop, AU, MU_SUN, DAY2SEC
    from math import pi, acos, cos, sin
    from scipy.linalg import norm

    legs = len(x) / 4
    common_mu = seq[0].mu_central_body

    # 1 -  we 'decode' the chromosome recording the various times of flight
    # (days) in the list T
    T = x[3::4]

    # 2 - We compute the epochs and ephemerides of the planetary encounters
    t_P = list([None] * (legs + 1))
    r_P = list([None] * (legs + 1))
    v_P = list([None] * (legs + 1))

    for i, planet in enumerate(seq):
        t_P[i] = epoch(start_mjd2000 + sum(T[:i]))
        r_P[i], v_P[i] = planet.eph(t_P[i])
        plot_planet(planet, t0=t_P[i], color=(
            0.8, 0.6, 0.8), legend=True, units = units, ax=axis)

    v_end_l = [a + b for a, b in zip(v_P[0], vinf_in)]
    # 4 - And we iterate on the legs
    for i in range(0, legs):
        # Fly-by
        v_out = fb_prop(v_end_l,
                        v_P[i],
                        x[1 + 4 * i] * seq[i].radius,
                        x[4 * i],
                        seq[i].mu_self)
        # s/c propagation before the DSM
        r, v = propagate_lagrangian(
            r_P[i], v_out, x[4 * i + 2] * T[i] * DAY2SEC, common_mu)
        plot_kepler(r_P[i], v_out, x[4 * i + 2] * T[i] * DAY2SEC,
                    common_mu, N=500, color='b', legend=False, units=units, ax=axis)
        # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
        dt = (1 - x[4 * i + 2]) * T[i] * DAY2SEC
        l = lambert_problem(r, r_P[i + 1], dt, common_mu, False, False)
        plot_lambert(
            l, sol=0, color='r', legend=False, units=units, N=500, ax=axis)
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]

# Plot of the trajectory of an mga_part problem


def _mga_part_plot(self, x):
    """
    Plots the trajectory represented by the decision vector x

    Example::

      prob.plot(x)
    """
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    axis = fig.gca(projection='3d', aspect='equal')

    # Plots the central 'planet'star
    axis.scatter(0, 0, 0, color='y')

    JR = 71492000.0
    seq = self.get_sequence()
    start_mjd2000 = self.t0.mjd2000
    _part_plot(x, JR, ax, seq, start_mjd2000, self.vinf_in)
    return ax
mga_part.plot = _mga_part_plot


def _mga_incipit_plot(self, x, plot_leg_0=False):
    """
    Plots the trajectory represented by the decision vector x

    Example::

      prob.plot(x)
    """
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from PyKEP.orbit_plots import plot_planet, plot_lambert, plot_kepler
    from PyKEP import epoch, propagate_lagrangian, lambert_problem, fb_prop, AU, MU_SUN, DAY2SEC
    from math import pi, acos, cos, sin
    from scipy.linalg import norm

    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d', aspect='equal')
    ax.scatter(0, 0, 0, color='y')

    JR = 71492000.0
    seq = self.get_sequence()
    common_mu = seq[0].mu_central_body
    r_P, v_P = seq[0].eph(epoch(x[0] + x[3]))

    # 3 - We start with the first leg: a lambert arc
    theta = 2 * pi * x[1]
    phi = acos(2 * x[2] - 1) - pi / 2
    # phi close to zero is in the moon orbit plane injection
    r = [cos(phi) * sin(theta), cos(phi) * cos(theta), sin(phi)]
    r = [JR * 1000 * d for d in r]

    l = lambert_problem(r, r_P, x[3] * DAY2SEC, common_mu, False, False)
    if (plot_leg_0):
        plot_lambert(ax, l, sol=0, color='k', legend=False, units=JR, N=500)

    # Lambert arc to reach seq[1]
    v_end_l = l.get_v2()[0]
    vinf_in = [a - b for a, b in zip(v_end_l, v_P)]
    _part_plot(x[4:], JR, ax, seq, x[0] + x[3], vinf_in)

    return ax
mga_incipit.plot = _mga_incipit_plot

# Plot of the trajectory for an mga_1dsm problem


def _mga_1dsm_tof_plot(self, x):
    """
    Plots the trajectory represented by the decision vector x
    """
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from PyKEP.orbit_plots import plot_planet, plot_lambert, plot_kepler
    from PyKEP import epoch, propagate_lagrangian, lambert_problem, fb_prop, AU, MU_SUN, DAY2SEC
    from math import pi, acos, cos, sin
    from scipy.linalg import norm

    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    axis = fig.gca(projection='3d')
    axis.scatter(0, 0, 0, color='y')

    seq = self.get_sequence()

    # 2 - We plot the first leg
    r_P0, v_P0 = seq[0].eph(epoch(x[0]))
    plot_planet(seq[0], t0=epoch(x[0]), color=(
        0.8, 0.6, 0.8), legend=True, units = AU, ax=axis)
    r_P1, v_P1 = seq[1].eph(epoch(x[0] + x[5]))
    theta = 2 * pi * x[1]
    phi = acos(2 * x[2] - 1) - pi / 2

    Vinfx = x[3] * cos(phi) * cos(theta)
    Vinfy = x[3] * cos(phi) * sin(theta)
    Vinfz = x[3] * sin(phi)

    v0 = [a + b for a, b in zip(v_P0, [Vinfx, Vinfy, Vinfz])]
    r, v = propagate_lagrangian(
        r_P0, v0, x[4] * x[5] * DAY2SEC, seq[0].mu_central_body)
    plot_kepler(
        r_P0,
        v0,
        x[4] *
        x[5] *
        DAY2SEC,
        seq[0].mu_central_body,
        N=100,
        color='b',
        legend=False,
        units=AU,
        ax=axis)

    # Lambert arc to reach seq[1]
    dt = (1 - x[4]) * x[5] * DAY2SEC
    l = lambert_problem(r, r_P1, dt, seq[0].mu_central_body)
    plot_lambert(l, sol=0, color='r', legend=False, units=AU, ax=axis)
    v_end_l = l.get_v2()[0]
    vinf_in = [a - b for a, b in zip(v_end_l, v_P1)]
    _part_plot(x[6:], AU, axis, seq[1:], x[0] + x[5], vinf_in)
    return axis
mga_1dsm_tof.plot = _mga_1dsm_tof_plot

del planet_js, planet_ss
