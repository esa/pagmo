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

import unittest as _ut


class _serialization_test(_ut.TestCase):

    def test_pickle(self):
        from PyGMO import archipelago, island_list, problem_list, algorithm_list, problem, algorithm
        import pickle
        from copy import deepcopy
        # We remove some problems that cannot be constructed without external
        # txt data files
        prob_list = deepcopy(problem_list)
        prob_list.remove(problem.cec2013)
        print('')
        for isl in island_list:
            for prob in prob_list:
                for algo in algorithm_list:
                    a = archipelago()
                    a.push_back(isl(algo(), prob(), 20))
                    a.push_back(isl(algo(), prob(), 20))
                    pickle.loads(pickle.dumps(a))

# This class will stress the island and archipelago classes with highly
# concurrent simple evolutions.


class _island_torture_test(_ut.TestCase):

    def __test_impl(self, isl_type, algo, prob):
        from PyGMO import archipelago, topology
        a = archipelago(topology=topology.ring())
        for i in range(0, 100):
            a.push_back(isl_type(algo, prob, 6))
        a.evolve(10)
        a.join()

    def test_local_island(self):
        from PyGMO import local_island, algorithm, problem
        isl_type = local_island
        algo_list = [algorithm.py_example(1), algorithm.de(5)]
        prob_list = [problem.py_example(), problem.dejong(1)]
        for algo in algo_list:
            for prob in prob_list:
                self.__test_impl(isl_type, algo, prob)

    def test_py_island(self):
        from PyGMO import py_island, algorithm, problem
        isl_type = py_island
        algo_list = [algorithm.py_example(1), algorithm.de(5)]
        prob_list = [problem.py_example(), problem.dejong(1)]
        for algo in algo_list:
            for prob in prob_list:
                self.__test_impl(isl_type, algo, prob)

    def test_ipy_island(self):
        from PyGMO import ipy_island, algorithm, problem
        try:
            from IPython.kernel import client
            mec = client.MultiEngineClient()
            if len(mec) == 0:
                raise RuntimeError()
        except ImportError as ie:
            return
        except BaseException as e:
            print(
                '\nThere is a problem with parallel IPython setup. The error message is:')
            print(e)
            print('Tests for ipy_island will not be run.')
            return
        isl_type = ipy_island
        algo_list = [algorithm.py_example(1), algorithm.de(5)]
        prob_list = [problem.py_example(), problem.dejong(1)]
        for algo in algo_list:
            for prob in prob_list:
                self.__test_impl(isl_type, algo, prob)


def run_serialization_test_suite():
    """Run the serialization test suite."""
    from PyGMO import test
    suite = _ut.TestLoader().loadTestsFromTestCase(_serialization_test)
    _ut.TextTestRunner(verbosity=2).run(suite)


def run_island_torture_test_suite():
    """Run the island torture test suite."""
    from PyGMO import test
    suite = _ut.TestLoader().loadTestsFromTestCase(_island_torture_test)
    _ut.TextTestRunner(verbosity=2).run(suite)


def run_full_test_suite():
    """Run the complete test suite for PyGMO."""
    from PyGMO import test
    from PyGMO.test._hypervolume_tests import get_hv_suite
    from PyGMO.test._topology_tests import get_topology_test_suite
    from PyGMO.test._archipelago_tests import get_archipelago_test_suite
    suite = _ut.TestLoader().loadTestsFromModule(test)

    # Add external suites explicitly
    suite.addTests(get_hv_suite())
    suite.addTests(get_topology_test_suite())
    suite.addTests(get_archipelago_test_suite())

    _ut.TextTestRunner(verbosity=2).run(suite)


def run_hv_test_suite():
    """Run the hypervolume test suite."""
    from PyGMO.test._hypervolume_tests import get_hv_suite
    _ut.TextTestRunner(verbosity=2).run(get_hv_suite())


def run_topology_test_suite():
    """Run the topology test suite."""
    from PyGMO.test._topology_tests import get_topology_test_suite
    _ut.TextTestRunner(verbosity=2).run(get_topology_test_suite())

def run_archipelago_test_suite():
    """Run the archipelago test suite."""
    from PyGMO.test._archipelago_tests import get_archipelago_test_suite
    _ut.TextTestRunner(verbosity=2).run(get_archipelago_test_suite())
