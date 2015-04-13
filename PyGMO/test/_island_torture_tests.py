# This class will stress the island and archipelago classes with highly
# concurrent simple evolutions.

import unittest
from PyGMO import archipelago, topology, local_island, algorithm, problem, py_island, ipy_island


class IslandTortureTests(unittest.TestCase):

    def __test_impl(self, isl_type, algo, prob):
        a = archipelago(topology=topology.ring())
        for i in range(0, 100):
            a.push_back(isl_type(algo, prob, 6))
        a.evolve(10)
        a.join()

    def test_local_island(self):
        isl_type = local_island
        algo_list = [algorithm.py_example(1), algorithm.de(5)]
        prob_list = [problem.py_example(), problem.dejong(1)]
        for algo in algo_list:
            for prob in prob_list:
                self.__test_impl(isl_type, algo, prob)

    def test_py_island(self):
        isl_type = py_island
        algo_list = [algorithm.py_example(1), algorithm.de(5)]
        prob_list = [problem.py_example(), problem.dejong(1)]
        for algo in algo_list:
            for prob in prob_list:
                self.__test_impl(isl_type, algo, prob)

    def test_ipy_island(self):
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

def get_island_torture_test_suite():
    suite = unittest.TestSuite()
    suite.addTests(unittest.makeSuite(IslandTortureTests))
    return suite
