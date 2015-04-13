from PyGMO import archipelago, island_list, problem_list, algorithm_list, problem
import pickle
from copy import deepcopy
import unittest

class SerializationTests(unittest.TestCase):

    def test_pickle(self):

        # We remove some problems that cannot be constructed without external
        # txt data files
        prob_list = deepcopy(problem_list)
        prob_list.remove(problem.cec2013)

        # Remove trajectory problems if PyKEP is not installed
        try:
            import PyKEP
        except ImportError:
            prob_list.remove(problem.py_pl2pl)

        print('')
        for isl in island_list:
            for prob in prob_list:
                for algo in algorithm_list:
                    a = archipelago()
                    a.push_back(isl(algo(), prob(), 20))
                    a.push_back(isl(algo(), prob(), 20))
                    pickle.loads(pickle.dumps(a))



def get_serialization_test_suite():
    suite = unittest.TestSuite()
    suite.addTests(unittest.makeSuite(SerializationTests))
    return suite
