import platform
import unittest
class IslandTortureTests(unittest.TestCase):

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

    @unittest.skipIf(platform.system() == "Windows", "The py_island test cannot be run on Windows")
    def test_py_island(self):
        # This test creates an infinite import loop on Windows when creating a new process
        # with Python multiprocessing module.
        # Normal usage of py_island on Windows is possible by using the `if __name__ == "__main__":` guard
        # Reference: http://stackoverflow.com/questions/20222534/python-multiprocessing-on-windows-if-name-main

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


def get_island_torture_test_suite():
    suite = unittest.TestSuite()
    suite.addTests(unittest.makeSuite(IslandTortureTests))
    return suite
