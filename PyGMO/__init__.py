# -*- coding: iso-8859-1 -*-
# Copyright (C) 2004-2009 The PaGMO development team,
# Advanced Concepts Team (ACT), European Space Agency (ESA)
# http://apps.sourceforge.net/mediawiki/pagmo
# http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers
# http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits
# act@esa.int
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

# Fill up the __extensions__ variable with all detected extensions
__extensions__ = {
    'nlopt': False,
    'gsl': False,
    'snopt': False,
    'ipopt': False,
    'gtop': False,
    'scipy': False,
    'networkx': False,
    'vpython': False,
    'pykep': False}
try:
    from scipy import __version__ as __scipy_version__
    __extensions__['scipy'] = True
except ImportError:
    pass
try:
    from networkx.version import version_info as __networkx_version__
    __extensions__['networkx'] = True
except ImportError:
    pass
try:
    from visual import version as __visual_version__
    __extensions__['vpython'] = True
except ImportError:
    pass
try:
    from PyKEP import __version__ as __pykep_version__
    __extensions__['pykep'] = True
except ImportError:
    pass

# NOTE: we may now test for extensions scipy, networkx, vpython and pykep, but not algos or problems
from PyGMO import core, algorithm, migration, problem, topology, test, util
from PyGMO.util import hypervolume, hv_algorithm

if "nlopt" in str(algorithm._get_algorithm_list()):
    __extensions__['nlopt'] = True
if "gsl" in str(algorithm._get_algorithm_list()):
    __extensions__['gsl'] = True
if "snopt" in str(algorithm._get_algorithm_list()):
    __extensions__['snopt'] = True
if "ipopt" in str(algorithm._get_algorithm_list()):
    __extensions__['ipopt'] = True

# NOTE: All extensions are now detected

__doc__ = 'PyGMO is a pretty cool guy. It kills aliens and does not afraid of anything...'
__all__ = ['core', 'algorithm', 'migration', 'problem', 'topology', 'test', 'util']
__version__ = '1.1.7dev'

# For convenience, bring all core classes into the root namespace when
# importing *.
from PyGMO.core import *
__all__ += [name for name in dir(core) if not name.startswith('_')]

problem_list = problem._get_problem_list()
algorithm_list = algorithm._get_algorithm_list()
island_list = core._get_island_list()


def run_test(n_trials=200, pop_size=20, n_gen=500):
    """
    This function runs some tests on the algorthm. Use it to verify the correct installation
    of PyGMO.

    USAGE: PyGMO.run_test(n_trials=200, pop_size = 20, n_gen = 500)

    * n_trials: each algorithm will be called n_trials times on the same problem to then evaluate best, mean and std
    * pop_size: this determines the population size
    * n_gen: this regulates the maximim number of function evaluation

    """
    from PyGMO import problem, algorithm, island
    from numpy import mean, std
    number_of_trials = n_trials
    number_of_individuals = pop_size
    number_of_generations = n_gen

    prob_list = [
        problem.schwefel(
            dim=10), problem.rastrigin(
            dim=10), problem.rosenbrock(
            dim=10), problem.ackley(
            dim=10), problem.griewank(
            dim=10), problem.levy5(10)]
    if __extensions__['gtop']:
        prob_list.append(problem.cassini_1())
        prob_list.append(problem.gtoc_1())
        prob_list.append(problem.cassini_2())
        prob_list.append(problem.messenger_full())

    algo_list = [
        algorithm.pso(
            gen=number_of_generations),
        algorithm.mde_pbx(
            gen=number_of_generations,
            xtol=1e-30,
            ftol=1e-30),
        algorithm.de(
            gen=number_of_generations,
            xtol=1e-30,
            ftol=1e-30),
        algorithm.jde(
            gen=number_of_generations,
            memory=False,
            xtol=1e-30,
            ftol=1e-30),
        algorithm.de_1220(
            gen=number_of_generations,
            memory=False,
            xtol=1e-30,
            ftol=1e-30),
        algorithm.sa_corana(
            iter=number_of_generations *
            number_of_individuals,
            Ts=1,
            Tf=0.01),
        algorithm.ihs(
            iter=number_of_generations *
            number_of_individuals),
        algorithm.sga(
            gen=number_of_generations),
        algorithm.cmaes(
            gen=number_of_generations,
            xtol=1e-30,
            ftol=1e-30,
            memory=False),
        algorithm.bee_colony(
            gen=number_of_generations /
            2)]
    print('\nTrials: ' + str(n_trials) + ' - Population size: ' +
          str(pop_size) + ' - Generations: ' + str(n_gen))
    for prob in prob_list:
        print('\nTesting problem: ' + prob.get_name() +
              ', Dimension: ' + str(prob.dimension))
        print('With Population Size: ' + str(pop_size))
        for algo in algo_list:
            print(' ' + str(algo))
            best = []
            best_x = []
            for i in range(0, number_of_trials):
                isl = island(algo, prob, number_of_individuals)
                isl.evolve(1)
                isl.join()
                best.append(isl.population.champion.f)
                best_x.append(isl.population.champion.x)
            print(' Best:\t' + str(min(best)[0]))
            print(' Mean:\t' + str(mean(best)))
            print(' Std:\t' + str(std(best)))

if __extensions__['scipy']:
    class race2algos:

        """
        This class uses the concept of racing to compare two algorithms
        on a probem. It runs repeatedly both algorithms on equal
        starting populations up to when it finds a statistical difference between
        the obtained samples. The difference is detected using Wilcoxon
        ranksum test. The algorithms are tested on populations of equal size.

        """

        def __init__(
                self,
                algo1,
                algo2,
                prob,
                pop_size=20,
                min_trials=20,
                p=0.05,
                max_runs=200):
            """
            Upon construction of the class object the race is initialized and launched.

            USAGE: r = PyGMO.race2algos(algo1, algo2, prob, pop_size=20, min_trials=20, p = 0.05, max_runs=200):

            * algo1: first algorithm in the race
            * algo2: second algorithm in the race
            * prob: problem (i.e. the "track" the algos are racing upon)
            * pop_size: population size of the island where the algos will perform evolution
            * min_trials: minimum number of runs to compare the algorithms
            * p: confidence level
            * max_runs: maximum number of races ....
            """
            from random import randint
            from copy import deepcopy
            from sys import stdout
            self.algo1 = algo1
            self.algo2 = algo2
            self.prob = prob
            self.res1 = []
            self.res2 = []
            self.pop_size = pop_size
            self.p = 0
            self.z = 0
            self.p_req = p
            print("Racing the algorithms ...")

            for i in range(max_runs):
                stdout.write("\rRuns: %i" % i)
                stdout.flush()
                # We reset the random number generators of the algorithm
                algo1.reset_rngs(randint(0, 9999999))
                algo2.reset_rngs(randint(0, 9999999))
                # We create an island with 20 individuals. This also initalizes
                # its population at random within the box bounds
                isl1 = island(algo1, prob, self.pop_size)
                # We copy the island and change its algo. Like this we make sure the two algorithms
                # will evolve the same inital population (good practice)
                isl2 = deepcopy(isl1)
                isl2.algorithm = algo2
                # We start the evolution (in parallel as we do not call the
                # method join())
                isl1.evolve(1)
                isl2.evolve(1)
                # Here join is called implicitly as we try to access one of the
                # islands during evolution
                self.res1.append(isl1.population.champion.f[0])
                self.res2.append(isl2.population.champion.f[0])
                # We check that the race is over (only after having accumulated
                # at least min_trials samples)
                if (i > min_trials):
                    if (self.are_different(self.res1, self.res2)):
                        break

        def are_different(self, data1, data2):
            from scipy.stats import wilcoxon
            self.z, self.p = wilcoxon(data1, data2)
            return (self.p < self.p_req)

        def plot(self):
            """
            Plots the result of the race

            USAGE: r.plot()
            """
            import matplotlib.pyplot as pl
            pl.subplot(1, 2, 1)
            pl.plot(sorted(self.res1), label="1." + self.algo1.get_name())
            pl.plot(sorted(self.res2), label="2." + self.algo2.get_name())
            pl.title(
                self.prob.get_name() + " dim: " + str(self.prob.dimension))
            pl.xlabel('rank')
            pl.legend()

            pl.subplot(1, 2, 2)
            pl.boxplot([self.res1, self.res2])
            pl.ylabel('Obj.Fun.')
            pl.title("Wilcoxon Test, p: %2.2e" %
                     self.p + ", z: " + str(self.z))
            pl.show()


def example_1(n_trials=25, variant_adptv=1, memory=True):
    from PyGMO import problem, algorithm, island, archipelago
    from PyGMO.topology import fully_connected
    from numpy import mean, median
    results = list()
    prob = problem.messenger_full()
    de_variants = [11, 13, 15, 17]
    algos = [
        algorithm.jde(
            gen=50,
            variant=v,
            memory=memory,
            variant_adptv=variant_adptv) for v in de_variants]

    for trial in range(n_trials):
        archi = archipelago(topology=fully_connected())
        for algo in algos:
            archi.push_back(island(algo, prob, 25))
        print("Trial N: " + str(trial))
        archi.evolve(30)
        results.append(min([isl.population.champion.f[0] for isl in archi]))
    return (mean(results), median(results), min(results), max(results))


def example_2(
        algo=algorithm.de(1),
        prob=problem.rosenbrock(10),
        topo=topology.barabasi_albert(
            3,
            3),
        n_evolve=100,
        n_isl=1024,
        pop_size=20,
        color_code='rank'):
    from PyGMO import problem, algorithm, island, archipelago
    from matplotlib.pyplot import savefig, close
    archi = archipelago(algo, prob, n_isl, pop_size, topology=topo)
    print("Drawing Initial Condition .. ")
    pos = archi.draw(
        scale_by_degree=True, n_size=3, e_alpha=0.03, n_color=color_code)
    savefig('archi000', dpi=72)
    close()
    for i in range(1, n_evolve):
        archi.evolve(1)
        archi.join()
        print("Drawing" + str(i) + "-th evolution .. ")
        pos = archi.draw(
            layout=pos,
            scale_by_degree=True,
            n_size=3,
            e_alpha=0.03,
            n_color=color_code)
        savefig('archi%03d' % i, dpi=72)
        close()
