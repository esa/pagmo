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
        from PyGMO import archipelago, island_list, problem_list, algorithm_list, problem
        import pickle
        from copy import deepcopy
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

# This class will stress the island and archipelago classes with highly
# concurrent simple evolutions.

import platform
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

    @_ut.skipIf(platform.system() == "Windows", "The py_island test cannot be run on Windows")
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
    import sys
    from PyGMO import test
    suite = _ut.TestLoader().loadTestsFromModule(test)

    # no need to Add external suites explicitly
    #suite.addTests(get_hv_suite())
    #suite.addTests(get_topology_test_suite())
    #suite.addTests(get_archipelago_test_suite())

    successful = _ut.TextTestRunner(verbosity=2).run(suite).wasSuccessful()
    sys.exit(0 if successful else 1)


def run_hv_test_suite():
    """Run the hypervolume test suite."""
    from PyGMO import test
    _ut.TextTestRunner(verbosity=2).run(get_hv_suite())


def run_topology_test_suite():
    """Run the topology test suite."""
    from PyGMO import test
    _ut.TextTestRunner(verbosity=2).run(get_topology_test_suite())


def run_archipelago_test_suite():
    """Run the archipelago test suite."""
    from PyGMO import test
    _ut.TextTestRunner(verbosity=2).run(get_archipelago_test_suite())


#Adding Archipalego test

class ArchipelagoTests(_ut.TestCase):
    """Test archipelago properties"""
    def test_simple_probs(self):
        """ Testing whether migration history matches the expected output for each migration_direction and distribution_type """
        from PyGMO import topology, algorithm, problem, archipelago, distribution_type, migration_direction, island, population, migration
        for migr_dir in [migration_direction.source, migration_direction.destination]:
            for dist_type in [distribution_type.point_to_point, distribution_type.broadcast]:
                prob = problem.rosenbrock(10)
                alg = algorithm.jde(20)
                archi = archipelago(alg, prob, 3, 20, migration_direction=migr_dir, distribution_type=dist_type)
                top = topology.ring(3)
                top.set_weight(0, 1, 0.0)
                top.set_weight(0, 2, 0.0)
                top.set_weight(1, 2, 0.0)
                top.set_weight(1, 0, 1.0)
                top.set_weight(2, 0, 1.0)
                top.set_weight(2, 1, 1.0)
                archi.topology = top
                archi.evolve(200)
                migr_hist = archi.dump_migr_history()
                # Below: After 200 evaluations, there should be some migrants from 1->0, 2->0 and 2->1
                # There should be no migrants from 1->0, 2->0 and 2->1
                self.assertTrue("(1,0,1)" not in migr_hist)
                self.assertTrue("(1,0,2)" not in migr_hist)
                self.assertTrue("(1,1,2)" not in migr_hist)

    def do_test_migr_setup(self, pop_xs, out_pop_xs, top, n_evolves):
        """ Generic procedure for testing whether the state of populations in 'pop_xs',
        after performing 'n_evolve' migration steps is equal to the expected 'out_pop_xs', given topology 'top'. """
        from PyGMO import topology, algorithm, problem, archipelago, distribution_type, migration_direction, island, population, migration
        prob = problem.identity()
        alg = algorithm.null()
        pops = []
        for xs in pop_xs:
            pop = population(prob)
            for x in xs:
                pop.push_back(x)
            pops.append(pop)

        archi = archipelago(distribution_type=distribution_type.broadcast, migration_direction=migration_direction.destination)
        for pop in pops:
            archi.push_back(island(alg, pop, s_policy=migration.best_s_policy(), r_policy=migration.fair_r_policy()))
        archi.topology = top
        archi.evolve_batch(n_evolves, 1, False)
        out_xs = []
        for ii, isl in enumerate(archi, 1):
            out_xs.append(tuple(sorted([i.cur_f for i in isl.population])))
        out_xs = tuple(out_xs)
        self.assertEqual(out_xs, out_pop_xs)

    def test_one_way_ring_null_alg(self):
        from PyGMO import topology
        pop_xs = (
            ((1.0, ), (2.0, ), (3.0, )),
            ((2.0, ), (2.0, ), (2.0, )),
            ((3.0, ), (3.0, ), (3.0, )),
            ((4.0, ), (4.0, ), (4.0, )),
        )
        # Since we evolve the populations in index-order, the champion of pop1 ([1, ])
        # will travel along the ring until pop4
        out_pop_xs = (
            ((1.0, ), (2.0, ), (3.0, )),
            ((1.0, ), (2.0, ), (2.0, )),
            ((1.0, ), (3.0, ), (3.0, )),
            ((1.0, ), (4.0, ), (4.0, )),
        )
        top = topology.one_way_ring(4)
        top.set_weight(1.0)
        self.do_test_migr_setup(pop_xs, out_pop_xs, top, 1)

        # We set the probability to 0.0 between islands 0 and 1. During migration,
        # it is the the champion of pop2 - (2.0, ) which travels along the ring
        top.set_weight(0, 1, 0.0)
        out_pop_xs_2 = (
            ((1.0, ), (2.0, ), (3.0, )),
            ((2.0, ), (2.0, ), (2.0, )),
            ((2.0, ), (3.0, ), (3.0, )),
            ((2.0, ), (4.0, ), (4.0, )),
        )
        self.do_test_migr_setup(pop_xs, out_pop_xs_2, top, 1)

        # We set the probability to 0.0, not migration should happen
        top.set_weight(0.0)
        self.do_test_migr_setup(pop_xs, pop_xs, top, 1)

    def test_distribution_type(self):
        from PyGMO import topology, algorithm, problem, archipelago, distribution_type, migration_direction, island, population, migration
        """Testing whether the distribution_type property works"""

        # Ensure the default distribution type is point_to_point
        a = archipelago()
        self.assertEqual(a.distribution_type, distribution_type.point_to_point)

        for dist_type in [distribution_type.point_to_point, distribution_type.broadcast]:
            a = archipelago()
            a.distribution_type = dist_type

            self.assertEqual(a.distribution_type, dist_type)

def get_archipelago_test_suite():
    suite = _ut.TestSuite()
    suite.addTests(_ut.makeSuite(ArchipelagoTests))
    return suite

#adding hypervolume tests

class HVCtorTest(_ut.TestCase):

    def setUp(self):
        self.good_ps_2d = [[5, 2], [3, 4], [1, 5]]

    def test_dimension_size(self):
        from PyGMO import hypervolume, hv_algorithm
        # fixed points dimensions should be reasonable
        # points of f_dim = 1
        self.assertRaises(ValueError, hypervolume, [[1, ], [2, ]])
        # points of f_dim = 0
        self.assertRaises(ValueError, hypervolume, [[], [], ])
        # empty set of points
        self.assertRaises(ValueError, hypervolume, [[], ])
        # point dimensions are not equal
        self.assertRaises(ValueError, hypervolume, [[1, 2, 3], [3, 3]])
        self.assertRaises(TypeError, hypervolume, [])  # empty set of points #2

    def test_pop_ctor(self):
        from PyGMO import hv_algorithm, problem, population, hypervolume
        # constructs the hypervolume object from a population object, expects
        # to not raise any error
        prob = problem.zdt(2)
        pop = population(prob, 100)

        # construction from a population object
        hypervolume(pop)

        # setting verification flag to False
        hypervolume(pop, False)

        # setting verification flag to True
        hypervolume(pop, True)

    def test_hypervolume_ctor_points(self):	
        from PyGMO import hv_algorithm, problem, population, hypervolume
        # test various possibilities of construction
        # something other than a list or a string
        self.assertRaises(TypeError, hypervolume, "A string")
        # good point definition with extra arg (hypervolume(ps, "something
        # extra")
        self.assertRaises(
            TypeError, hypervolume, [[1, 2, 3], [2, 3, 4], [2, 3, 5]],
            "extra arg")
        self.assertRaises(
            TypeError, hypervolume, [[1, 2, 3], [2, 3, 4], "foo"])  # bad point
        # bad point value
        self.assertRaises(
            TypeError, hypervolume, [[1, 2, 3], [2, 3, 4], [2, 3, "bar"]])
        # skipping argument: hypervolume() should raise TypeError
        self.assertRaises(TypeError, hypervolume)

        self.assertRaises(
            TypeError, hypervolume, foo=self.good_ps_2d)  # bad kwarg
        self.assertRaises(
            TypeError, hypervolume, self.good_ps_2d, foo="bar")  # extra kwarg


class HVFlagsTest(_ut.TestCase):

    def test_gettersetter(self):
        from PyGMO import hv_algorithm, hypervolume
        hv = hypervolume([[1, 2, 3], [4, 5, 6]])
        self.assertTrue(hv.get_verify())
        self.assertTrue(hv.get_copy_points())
        hv.set_verify(False)
        hv.set_copy_points(False)
        self.assertTrue(hv.get_verify() is False)
        self.assertTrue(hv.get_copy_points() is False)


class HVComputeTest(_ut.TestCase):

    def setUp(self):
	from PyGMO import hypervolume, hv_algorithm
        self.hv2d = hypervolume([[3, 1], [2, 2], [1, 3]])

    def test_correct_out(self):
	from PyGMO import hypervolume, hv_algorithm
        # simple 3D test
        hv = hypervolume([[1, 1, 1], [2, 2, 2, ]])
        self.assertEqual(hv.compute(r=[3, 3, 3]), 8)

        # simple 2D test
        hv = hypervolume([[1, 2], [2, 1]])
        self.assertEqual(hv.compute(r=[3, 3]), 3)

        # point on the border of refpoint (2D)
        hv = hypervolume([[1, 2], [2, 1]])
        self.assertEqual(hv.compute([2, 2]), 0)

        # points on the border of refpoint (3D)
        hv = hypervolume([[1, 2, 1], [2, 1, 1]])
        self.assertEqual(hv.compute([2, 2, 2]), 0)

    def test4d_dominated(self):
        from PyGMO import hypervolume, hv_algorithm
        hv = hypervolume([[1, 1, 1, 1], [2, 2, 2, 2]])
        self.assertEqual(
            hv.compute(r=(3, 3, 3, 3), algorithm=hv_algorithm.hv4d()), 16)
        self.assertEqual(
            hv.compute(r=(3, 3, 3, 3), algorithm=hv_algorithm.fpl()), 16)

    def test4d_edge(self):
        from PyGMO import hypervolume, hv_algorithm
        hv = hypervolume([[1, 1, 1, 3], [2, 2, 2, 3]])
        self.assertEqual(
            hv.compute(r=(3, 3, 3, 3), algorithm=hv_algorithm.hv4d()), 0)
        self.assertEqual(
            hv.compute(r=(3, 3, 3, 3), algorithm=hv_algorithm.fpl()), 0)

    def test4d_duplicate(self):
        from PyGMO import hypervolume, hv_algorithm
        hv = hypervolume([[1, 1, 1, 1], [1, 1, 1, 1]])
        self.assertEqual(
            hv.compute(r=(2, 2, 2, 2), algorithm=hv_algorithm.hv4d()), 1)
        self.assertEqual(
            hv.compute(r=(2, 2, 2, 2), algorithm=hv_algorithm.fpl()), 1)

        # Duplicate and dominated
        hv = hypervolume([[1, 1, 1, 1], [1, 1, 1, 1], [0, 0, 0, 0]])
        self.assertEqual(
            hv.compute(r=(2, 2, 2, 2), algorithm=hv_algorithm.hv4d()), 16)
        self.assertEqual(
            hv.compute(r=(2, 2, 2, 2), algorithm=hv_algorithm.fpl()), 16)

    def test_tuple_ctor(self):
        # test that hypervolume can be computed using a tuple as well
	from PyGMO import hypervolume, hv_algorithm
        hv = hypervolume(((1, 1, 1), (2, 2, 2,)))
        self.assertEqual(hv.compute(r=(3, 3, 3)), 8)

    def test_casting_float(self): 
        # casting to float
        from PyGMO import hypervolume, hv_algorithm
        self.assertEqual(self.hv2d.compute(["4.0", "4"]), 6)

    def test_refpoint_not_dom(self):
        # refpoint must be at least weakly dominated by every point (assuming
        # minimization problem)
        from PyGMO import hypervolume, hv_algorithm
        hv = hypervolume([[1, 3], [2, 2], [3, 1]])
        # equal to some other point
        self.assertRaises(ValueError, hv.compute, [3, 1])
        # refpoint dominating some points
        self.assertRaises(ValueError, hv.compute, [1.5, 1.5])
        # refpoint dominating all points
        self.assertRaises(ValueError, hv.compute, [0, 0])

    def test_kwargs(self):
        # using kwarg 'r' correctly
        from PyGMO import hypervolume, hv_algorithm
        self.assertEqual(self.hv2d.compute(r=[3.5, 3.5]), 3.25)
        # using kwargs 'r', 'algorithm' correctly
        self.assertEqual(
            self.hv2d.compute(r=[3.5, 3.5], algorithm=hv_algorithm.hv2d()),
            3.25)

        # bad kwarg for reference point
        self.assertRaises(TypeError, self.hv2d.compute, refpoint=[4, 4])
        # we do not accept random kwargs
        self.assertRaises(
            TypeError, self.hv2d.compute, r=[4, 4], foo="Something extra")
        # we do not accept random kwargs (as above but with ref point as arg)
        self.assertRaises(
            TypeError, self.hv2d.compute, [4, 4], foo="Something extra")
        self.assertRaises(
            TypeError, self.hv2d.compute, [4, 4], hv_algorithm.hv2d(),
            foo="Something extra")   # we do not accept random kwargs

    def test_kwargs_hv_algo(self):
        # using kwargs correctly
        from PyGMO import hypervolume, hv_algorithm
        self.assertEqual(
            self.hv2d.compute(r=[4, 4], algorithm=hv_algorithm.hv2d()), 6)
        # using kwargs in reversed order
        self.assertEqual(
            self.hv2d.compute(algorithm=hv_algorithm.hv2d(), r=[4, 4]), 6)
        # arg + kwarg
        self.assertEqual(
            self.hv2d.compute([4, 4], algorithm=hv_algorithm.hv2d()), 6)
        self.assertEqual(
            self.hv2d.compute([4, 4], hv_algorithm.hv2d()), 6)  # arg + arg
        # should use nadir point
        self.assertRaises(
            TypeError, self.hv2d.compute, algorithm=hv_algorithm.hv2d())

    def test_bad_algo(self):
        from PyGMO import hypervolume, hv_algorithm
        # 3d method to 2d problem
        self.assertRaises(
            ValueError, self.hv2d.compute, [4, 4], hv_algorithm.hv3d())


class HVContributionsTest(_ut.TestCase):

    def assertContribs(self, S, R, ans):
        from PyGMO import hypervolume, hv_algorithm
        """
        This method is an assertion that given hypervolume problem constructed from S, with a reference point R
        Returns a valid answer to the "contributions" feature both for the contributions method and the explicit call
        for the exclusive hypervolume as well.
        """
        hv = hypervolume(S)
        self.assertEqual(hv.contributions(R), ans)
        self.assertEqual(tuple(hv.exclusive(i, R) for i in range(len(S))), ans)

    def test2d(self):
        """
        This test contains a front with 3 non dominated points,
        and many dominated points. Most of the dominated points
        lie on edges of the front, which makes their exclusive contribution
        equal to 0.
        """
        S = ((1, 6.5), (1, 6), (1, 5), (2, 5), (3, 5), (3, 3), (4, 6.5),
             (4.5, 4), (5, 3), (5, 1.5), (7, 1.5), (7, 3.5), )
        R = (7, 6.5, )
        ans = (0.0, 0.0, 1.0, 0.0, 0.0, 3.5, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, )

        self.assertContribs(S, R, ans)

        # Adding few extra points that share an edge with a reference point
        extra = ((7, 0.5), (7, 1.0), (7, 4.5), (0.0, 6.5), (5.5, 6.5), )
        S += extra
        ans += (0, ) * len(extra)

        self.assertContribs(S, R, ans)

        # Adding few duplicated points on the inside and on the edges
        extra = ((7, 0.5), (5.5, 6.5), (5, 5), (5, 5), (5, 5), )
        S += extra
        ans += (0, ) * len(extra)

        self.assertContribs(S, R, ans)

    def test2d_gradual(self):
        """
        Gradually adding duplicate points to the set, making sure the contribution change accordingly.
        """
        S = ((1, 1),)
        R = (2, 2)
        self.assertContribs(S, R, (1,))

        S += ((1, 1),)
        self.assertContribs(S, R, (0, 0))

        S += ((1, 1),)
        self.assertContribs(S, R, (0, 0, 0))

        S += ((0.5, 0.5),)
        self.assertContribs(S, R, (0, 0, 0, 1.25))

        S += ((0.5, 0.5),)
        self.assertContribs(S, R, (0,) * 5)

    def test3d(self):
        """
        This test contains a tricky front in 3D with some weakly dominated points on the "edges" of the bounding box.
        """

        # Non-tricky base problem
        S = ((-6, -1, -6), (-1, -3, -5), (-3, -4, -4),
             (-4, -2, -3), (-5, -5, -2), (-2, -6, -1),)
        R = (0, 0, 0)
        ans = (18, 2, 12, 1, 18, 2)

        self.assertContribs(S, R, ans)

        # Add some points that contribute nothing and do not alter other
        # contributions
        extra = ((-3, -1, -3), (-1, -1, -5), (-1, -2, -4),
                 (-1, -3, -4), (-7, -7, 0), (0, -5, -5), (-7, 0, -7))

        S += extra
        ans += (0,) * len(extra)
        self.assertContribs(S, R, ans)

    def test3d_gradual(self):
        """
        Gradually adding points, some of which are dominated or duplicates.
        Tests whether contributions and repeated exclusive method produce the same results.
        """
        S = ((3, 3, 3),)
        R = (5, 5, 5)
        self.assertContribs(S, R, (8,))

        # Decrease the contribution of first point. Second point is dominated.
        S += ((4, 4, 4),)
        self.assertContribs(S, R, (7, 0,))

        # Add duplicate point
        S += ((3, 3, 3),)
        self.assertContribs(S, R, (0, 0, 0))

        S += ((3, 3, 2),)
        self.assertContribs(S, R, (0, 0, 0, 4))

        S += ((3, 3, 1),)
        self.assertContribs(S, R, (0, 0, 0, 0, 4))

    def test3d_extreme(self):
        from PyGMO import hypervolume
        """
        Combine extreme points together.
        Mixing small and large contributions in a single front
        """

        # Reset the set S.
        # 3 duplicate points
        R = (0, 0, 0)
        S = ((-1, -1, -1),) * 3
        self.assertContribs(S, R, (0,) * 3)

        # Adding a point far away
        S += ((-1000,) * 3,)
        self.assertContribs(S, R, (0, 0, 0, 999999999))

        # Adding an even further point
        S += ((-10000,) * 3,)
        self.assertContribs(S, R, (0, 0, 0, 0, 999000000000))

        # Tiny box on top of a large one
        S = ((-1000.001, -0.001, -0.001), (-1000, -1000, -1000))
        R = (0, 0, 0)
        hv = hypervolume(S)
        ans = (0.000000001, 999999999.999)
        c = list(hv.contributions(R))
        # Round contribution to 9th decimal place as the double type is loosing
        # the exact accuracy
        c[0] = round(c[0], 9)
        self.assertEqual(tuple(c), ans)

    def test4d(self):
        """
        Gradually adding points.
        Tests whether contributions and repeated exclusive methods produce the same results.
        """
        S = ((1, 1, 1, 1), )
        R = (5, 5, 5, 5)
        self.assertContribs(S, R, (256,))

        S += ((4, 4, 4, 4),)
        self.assertContribs(S, R, (255, 0, ))

        S += ((3, 3, 3, 3),)
        self.assertContribs(S, R, (240, 0, 0))

        S += ((1, 1, 1, 1),)
        self.assertContribs(S, R, (0, ) * 4)

    def test5d(self):
        """
        Gradually adding points.
        Tests whether contributions and repeated exclusive methods produce the same results.
        """
        S = ((1, 1, 1, 1, 1), )
        R = (5, 5, 5, 5, 5)
        self.assertContribs(S, R, (1024,))

        S += ((4, 4, 4, 4, 4),)
        self.assertContribs(S, R, (1023, 0, ))

        S += ((3, 3, 3, 3, 3),)
        self.assertContribs(S, R, (992, 0, 0,))

        S += ((1, 1, 1, 1, 1),)
        self.assertContribs(S, R, (0,) * 4)


class HVLeastContribTest(_ut.TestCase):

    def setUp(self):
        from PyGMO import hypervolume, hv_algorithm
        self.r = [4, 4]
        self.hv2d_eq_0 = hypervolume([[3, 1], [2, 2], [1, 3]])  # LC in [0,1,2]
        self.hv2d_eq_1 = hypervolume([[2.5, 1], [2, 2], [1, 3]])  # LC = 1
        self.hv2d_0 = hypervolume([[3.5, 1], [2, 2], [1, 3]])
        self.hv2d_1 = hypervolume([[3, 1], [2.5, 2.5], [1, 3]])
        self.hv2d_2 = hypervolume([[3, 1], [2, 2], [1, 3.5]])

    def test_correct_out(self):
        from PyGMO import hypervolume, hv_algorithm
        self.assertTrue(
            self.hv2d_eq_0.least_contributor(r=self.r) in [0, 1, 2])

        self.assertEqual(self.hv2d_eq_1.least_contributor(r=self.r), 1)
        self.assertEqual(self.hv2d_0.least_contributor(r=self.r), 0)
        self.assertEqual(self.hv2d_1.least_contributor(r=self.r), 1)
        self.assertEqual(self.hv2d_2.least_contributor(r=self.r), 2)

    def test_kwargs(self):
        from PyGMO import hypervolume, hv_algorithm
        # using kwarg 'r' correctly
        self.assertEqual(self.hv2d_1.least_contributor(r=self.r), 1)
        # using kwarg 'r' and 'algorithm' correctly
        self.assertEqual(
            self.hv2d_1.least_contributor(
                r=self.r,
                algorithm=hv_algorithm.hv2d()),
            1)

        # bad kwarg for reference point
        self.assertRaises(
            TypeError, self.hv2d_0.least_contributor, refpoint=[4, 4])
        self.assertRaises(
            TypeError, self.hv2d_0.least_contributor, r=[
                4, 4], foo="Something extra")  # we do not accept random kwargs
        # we do not accept random kwargs (as above but with ref point as arg)
        self.assertRaises(
            TypeError, self.hv2d_0.least_contributor, [
                4, 4], foo="Something extra")
        self.assertRaises(
            TypeError, self.hv2d_0.least_contributor, [4, 4], hv_algorithm.
            hv2d(), foo="Something extra")   # we do not accept random kwargs

    def test_bad_algo(self):
        from PyGMO import hypervolume, hv_algorithm
        self.assertRaises(
            ValueError, self.hv2d_0.least_contributor, [
                4, 4], hv_algorithm.hv3d())  # 3d method to 2d problem


class HVExclusiveTest(_ut.TestCase):

    def setUp(self):
        from PyGMO import hypervolume, hv_algorithm
        self.r = [4, 4]
        # all are equal (take first -> idx = 0)
        self.hv2d = hypervolume([[3, 1], [2, 2], [1, 3]])
        # all are equal (take first -> idx = 0)
        self.hv2d_2 = hypervolume([[3.1, 1], [2, 2], [1, 3]])

    def test_correct_out(self):
        from PyGMO import hypervolume, hv_algorithm
        self.assertEqual(self.hv2d.exclusive(p_idx=0, r=self.r), 1)
        self.assertEqual(self.hv2d.exclusive(p_idx=1, r=self.r), 1)
        self.assertEqual(self.hv2d.exclusive(p_idx=2, r=self.r), 1)
        self.assertTrue(
            abs(self.hv2d_2.exclusive(p_idx=0, r=self.r) - 0.9) < 0.00000001)

    def test_kwargs(self):
        from PyGMO import hypervolume, hv_algorithm
        self.assertRaises(TypeError, self.hv2d.exclusive, 0)  # no refpoint
        # no refpoint
        self.assertRaises(TypeError, self.hv2d.exclusive, p_idx=0)
        self.assertRaises(
            TypeError,
            self.hv2d.exclusive,
            p_idx=0,
            hv_algorithm=hv_algorithm.wfg())  # no refpoint
        self.assertEqual(self.hv2d.exclusive(0, self.r), 1)  # using arg
        # using kwarg 'r' correctly
        self.assertEqual(self.hv2d.exclusive(0, r=self.r), 1)
        # using kwarg 'r' correctly
        self.assertEqual(self.hv2d.exclusive(p_idx=0, r=self.r), 1)

        self.assertRaises(
            TypeError,
            self.hv2d.exclusive,
            p_idx=0,
            algorithm=hv_algorithm.hv2d())  # no refpoint
        self.assertEqual(
            self.hv2d.exclusive(0, self.r, hv_algorithm.hv2d()), 1)  # all args
        self.assertEqual(
            self.hv2d.exclusive(
                0,
                self.r,
                algorithm=hv_algorithm.hv2d()),
            1)  # last kwarg
        self.assertEqual(self.hv2d.exclusive(
            p_idx=0, r=self.r, algorithm=hv_algorithm.hv2d()), 1)  # all kwargs
        self.assertEqual(
            self.hv2d.exclusive(
                algorithm=hv_algorithm.hv2d(),
                r=self.r,
                p_idx=0),
            1)  # all kwargs in reverse

        # bad kwarg for reference point
        self.assertRaises(TypeError, self.hv2d.exclusive, 0, refpoint=[4, 4])
        # we do not accept random kwargs
        self.assertRaises(
            TypeError, self.hv2d.exclusive, 0, r=[4, 4], foo="Something extra")
        # we do not accept random kwargs (as above but with ref point as arg)
        self.assertRaises(
            TypeError, self.hv2d.exclusive, 0, [4, 4], foo="Something extra")
        self.assertRaises(
            TypeError, self.hv2d.exclusive, 0, [4, 4], hv_algorithm.hv2d(),
            foo="Something extra")   # we do not accept random kwargs
        self.assertRaises(
            TypeError,
            self.hv2d.exclusive,
            p_idx=0,
            r=self.r,
            algorithm=hv_algorithm.hv2d(),
            foo="Something extra")
        self.assertRaises(TypeError, self.hv2d.exclusive,
                          r=self.r, algorithm=hv_algorithm.hv2d())  # no p_idx

    def test_p_idx(self):
        from PyGMO import hypervolume, hv_algorithm
        # negative idx
        self.assertRaises(TypeError, self.hv2d.exclusive, -1, self.r)
        self.assertRaises(
            ValueError, self.hv2d.exclusive, 100, self.r)  # large
        self.assertRaises(
            TypeError, self.hv2d.exclusive, "not an int", self.r)  # large
        self.assertRaises(TypeError, self.hv2d.exclusive)  # p_idx not provided

    def test_bad_algo(self):
	from PyGMO import hypervolume, hv_algorithm
        # 3d method to 2d problem
        self.assertRaises(
            ValueError, self.hv2d.exclusive, 0, [4, 4], hv_algorithm.hv3d())


class HVNadirPointTest(_ut.TestCase):

    def setUp(self):
	from PyGMO import hypervolume, hv_algorithm        
	self.hv2d = hypervolume([[3, 1], [2, 2], [1, 3]])

    def test_nadir_point(self):
	from PyGMO import hypervolume, hv_algorithm
        # default nadir point
        self.assertEqual(tuple(self.hv2d.get_nadir_point()), (3, 3))
        # custom nadir point
        self.assertEqual(tuple(self.hv2d.get_nadir_point(5.0)), (8, 8))
        # nadir point with eps=0.0
        self.assertEqual(tuple(self.hv2d.get_nadir_point(0.0)), (3, 3))
        # nadir point with eps=-0.0 is ok
        self.assertEqual(tuple(self.hv2d.get_nadir_point(-0.0)), (3, 3))
        # custom nadir point with 'eps' kwarg
        self.assertEqual(tuple(self.hv2d.get_nadir_point(eps=5.0)), (8, 8))
        # nadir point with negative eps
        self.assertRaises(ValueError, self.hv2d.get_nadir_point, -0.0000001)
        self.assertRaises(
            TypeError, self.hv2d.get_nadir_point, eps="foo")  # bad kwarg
        self.assertRaises(
            TypeError, self.hv2d.get_nadir_point, "foo")  # bad arg
        # bad kwarg name
        self.assertRaises(TypeError, self.hv2d.get_nadir_point, epsilon=1.0)


class HVAlgorithms(_ut.TestCase):

    def setUp(self):
	from PyGMO import hypervolume, hv_algorithm
        self.ps4d_0 = [[1., ] * 4, ]
        self.ps3d_0 = [[1., ] * 3, ]
        self.r4d_0 = [2, ] * 4
        self.r3d_0 = [2, ] * 3

    def test_wfg(self):
	from PyGMO import hypervolume, hv_algorithm
        # stop_dimension = 0
        self.assertRaises(ValueError, hv_algorithm.wfg, stop_dimension=0)
        # stop_dimension = 1
        self.assertRaises(ValueError, hv_algorithm.wfg, stop_dimension=1)

    def test_hv4d(self):
	from PyGMO import hypervolume, hv_algorithm
        hv3d = hypervolume(self.ps3d_0)
        hv4d = hypervolume(self.ps4d_0)

        self.assertEqual(
            1.0, hv4d.compute(self.r4d_0, algorithm=hv_algorithm.hv4d()))
        self.assertRaises(
            ValueError,
            hv3d.compute,
            r=self.r3d_0,
            algorithm=hv_algorithm.hv4d())
        self.assertRaises(
            ValueError,
            hv3d.compute,
            r=self.r4d_0,
            algorithm=hv_algorithm.hv4d())
        self.assertRaises(
            ValueError,
            hv4d.compute,
            r=self.r3d_0,
            algorithm=hv_algorithm.hv4d())


def get_hv_suite():
    suite = _ut.TestSuite()
    suite.addTests(_ut.makeSuite(HVCtorTest))
    suite.addTests(_ut.makeSuite(HVComputeTest))
    suite.addTests(_ut.makeSuite(HVLeastContribTest))
    suite.addTests(_ut.makeSuite(HVExclusiveTest))
    suite.addTests(_ut.makeSuite(HVNadirPointTest))
    suite.addTests(_ut.makeSuite(HVAlgorithms))
    suite.addTests(_ut.makeSuite(HVFlagsTest))
    suite.addTests(_ut.makeSuite(HVContributionsTest))
    return suite

#adding topology tests

class TopologyTests(_ut.TestCase):
    """ Tests for topology object and migration probabability """

    def test_overloaded_add_edge(self):
        from PyGMO import topology, algorithm, problem, archipelago, distribution_type, migration_direction, island, population, migration
        """ Tests the add_edge(n, m, weight) setter """

        t = topology.custom()
        t.push_back()
        t.push_back()
        t.push_back()
        t.add_edge(0, 1, 0.1)
        t.add_edge(2, 0)  # 1.0 by default
        t.add_edge(0, 2, 0.3)
        self.assertEqual(t.get_weight(0, 1), 0.1)
        self.assertEqual(t.get_weight(2, 0), 1.0)
        self.assertEqual(t.get_weight(0, 2), 0.3)

    def test_incorrect_probability(self):
        from PyGMO import topology, algorithm, problem, archipelago, distribution_type, migration_direction, island, population, migration
        """ Tests the setter for probabilities outside [0, 1] """

        top = topology.custom()
        top.push_back()
        top.push_back()

        # Add edges with valid probability
        top.add_edge(0, 1, 0.1)
        top.add_edge(1, 0, 0.1)
        top.set_weight(0.1)
        top.set_weight(0, 0.1)
        top.set_weight(0, 1, 0.1)

        self.assertRaises(ValueError, top.set_weight, -0.1)
        self.assertRaises(ValueError, top.set_weight, 0, -0.1)
        self.assertRaises(ValueError, top.set_weight, 0, 1, -0.1)

        self.assertRaises(ValueError, top.set_weight, 1.1)
        self.assertRaises(ValueError, top.set_weight, 0, 1.1)
        self.assertRaises(ValueError, top.set_weight, 0, 1, 1.1)

        # Add new vertex and try to add_edge with incorrect probability
        top.push_back()
        self.assertRaises(ValueError, top.add_edge, 1, 2, -0.1)
        self.assertRaises(ValueError, top.add_edge, 2, 1, 1.1)

    def test_basic_probs(self):
        from PyGMO import topology, algorithm, problem, archipelago, distribution_type, migration_direction, island, population, migration
        """ Tests the basic getter/setter for custom migration probability """
        t = topology.custom()
        t.push_back()
        t.push_back()
        t.push_back()
        t.add_edge(0, 1)
        t.add_edge(2, 0)
        t.add_edge(0, 2)
        t.set_weight(2, 0, 0.1)
        t.set_weight(0, 2, 0.2)
        self.assertEqual(t.get_weight(0, 1), 1.0)  # This one should be set by default
        self.assertEqual(t.get_weight(2, 0), 0.1)
        self.assertEqual(t.get_weight(0, 2), 0.2)

    def test_ring_topo(self):
        from PyGMO import topology, algorithm, problem, archipelago, distribution_type, migration_direction, island, population, migration
        """ Tests the getter/setter for the ring topology """

        t = topology.ring()
        t.push_back()
        t.push_back()

        # By default the weights are 1.0
        self.assertEqual(t.get_weight(0, 1), 1.0)
        self.assertEqual(t.get_weight(1, 0), 1.0)

        # Once you alter one weight, the remaining are still 1.0
        t.set_weight(0, 1, 0.1)
        self.assertEqual(t.get_weight(0, 1), 0.1)
        self.assertEqual(t.get_weight(1, 0), 1.0)

        # Pushing new node should not change anything to the previous ones
        t.push_back()
        self.assertEqual(t.get_weight(0, 1), 0.1)
        self.assertEqual(t.get_weight(1, 0), 1.0)
        self.assertEqual(t.get_weight(0, 2), 1.0)
        self.assertEqual(t.get_weight(2, 0), 1.0)
        self.assertEqual(t.get_weight(1, 2), 1.0)
        self.assertEqual(t.get_weight(2, 1), 1.0)

        # Test the all-custom weights for ring topology
        t.set_weight(0, 1, 0.1)
        t.set_weight(1, 0, 0.2)
        t.set_weight(0, 2, 0.3)
        t.set_weight(2, 0, 0.4)
        t.set_weight(1, 2, 0.5)
        t.set_weight(2, 1, 0.6)
        self.assertEqual(t.get_weight(0, 1), 0.1)
        self.assertEqual(t.get_weight(1, 0), 0.2)
        self.assertEqual(t.get_weight(0, 2), 0.3)
        self.assertEqual(t.get_weight(2, 0), 0.4)
        self.assertEqual(t.get_weight(1, 2), 0.5)
        self.assertEqual(t.get_weight(2, 1), 0.6)

    def test_mixed_migr_representation(self):
        from PyGMO import topology, algorithm, problem, archipelago, distribution_type, migration_direction, island, population, migration
        """ Testing the mixing of different edge weight setters """

        t = topology.custom()
        t.push_back()
        t.push_back()
        t.push_back()
        t.push_back()
        t.add_edge(0, 1)
        t.add_edge(1, 0)
        t.add_edge(1, 2)
        self.assertEqual(t.get_weight(0, 1), 1.0)
        self.assertEqual(t.get_weight(1, 0), 1.0)
        self.assertEqual(t.get_weight(1, 2), 1.0)
        t.set_weight(1, 0.5)  # Set weight for each out-going edge of vertex 1
        self.assertEqual(t.get_weight(0, 1), 1.0)
        self.assertEqual(t.get_weight(1, 0), 0.5)
        self.assertEqual(t.get_weight(1, 2), 0.5)
        t.add_edge(1, 3)  # Weight of a new out-going edge of vertex 1 should still be 1.0 by default
        self.assertEqual(t.get_weight(1, 0), 0.5)
        self.assertEqual(t.get_weight(1, 2), 0.5)
        self.assertEqual(t.get_weight(1, 3), 1.0)

        t.set_weight(0.3)  # Set the weight to all edges
        self.assertEqual(t.get_weight(0, 1), 0.3)
        self.assertEqual(t.get_weight(1, 0), 0.3)
        self.assertEqual(t.get_weight(1, 2), 0.3)
        self.assertEqual(t.get_weight(1, 3), 0.3)
        t.add_edge(2, 3)  # Weight of a new edge should still be 1.0 by default

        self.assertEqual(t.get_weight(0, 1), 0.3)
        self.assertEqual(t.get_weight(1, 0), 0.3)
        self.assertEqual(t.get_weight(1, 2), 0.3)
        self.assertEqual(t.get_weight(1, 3), 0.3)
        self.assertEqual(t.get_weight(2, 3), 1.0)

    def test_chain_break(self):
        from PyGMO import topology, algorithm, problem, archipelago, distribution_type, migration_direction, island, population, migration
        """ Testing 'chain breaking' after push_back in one-way-ring """

        t = topology.ring(2)
        t.set_weight(0, 1, 0.1)
        t.set_weight(1, 0, 0.2)
        t.push_back()
        # In two-way ring, the weights between the first and second island should be the same after push_back
        self.assertEqual(t.get_weight(0, 1), 0.1)
        self.assertEqual(t.get_weight(1, 0), 0.2)
        self.assertEqual(t.get_weight(1, 2), 1.0)
        self.assertEqual(t.get_weight(2, 1), 1.0)
        self.assertEqual(t.get_weight(2, 0), 1.0)
        self.assertEqual(t.get_weight(0, 2), 1.0)

        # In one-way-ring, the edge 1->0 will be broken
        # and replaced by 1->2 and 2->0 with weights 1.0 each
        t = topology.one_way_ring(2)
        t.set_weight(0, 1, 0.1)
        t.set_weight(1, 0, 0.2)
        t.push_back()
        self.assertEqual(t.get_weight(0, 1), 0.1)
        self.assertEqual(t.get_weight(1, 2), 1.0)
        self.assertEqual(t.get_weight(2, 0), 1.0)

    def test_topology_serialize(self):
        from PyGMO import topology, algorithm, problem, archipelago, distribution_type, migration_direction, island, population, migration
        """ Testing whether the weights are retained after serialization of an archipelago """

        prob = problem.rosenbrock(10)
        alg = algorithm.jde(20)
        archi = archipelago(alg, prob, 4, 20)
        top = topology.ring(4)
        top.set_weight(0, 1, 0.01)
        top.set_weight(0, 3, 0.03)
        top.set_weight(1, 2, 0.12)
        top.set_weight(2, 1, 0.21)
        archi.topology = top
        archi.evolve(5)
        import pickle
        pickle.loads(pickle.dumps(archi))
        self.assertEqual(archi.topology.get_weight(0, 1), 0.01)
        self.assertEqual(archi.topology.get_weight(0, 3), 0.03)
        self.assertEqual(archi.topology.get_weight(1, 2), 0.12)
        self.assertEqual(archi.topology.get_weight(2, 1), 0.21)

    def test_incorrect_weights(self):
        from PyGMO import topology, algorithm, problem, archipelago, distribution_type, migration_direction, island, population, migration
        """ Testing ValueError for getter/setter """

        top = topology.ring(4)
        self.assertRaises(ValueError, top.get_weight, 0, 2)
        self.assertRaises(ValueError, top.set_weight, 0, 2, 0.5)

    def test_average_path_length(self):
	from PyGMO import topology
        """ Testing the average path length algorithm """

        top = topology.ring(2)
        self.assertEqual(top.get_average_shortest_path_length(), 1.0)
        top = topology.ring(3)
        self.assertEqual(top.get_average_shortest_path_length(), 1.0)
        top = topology.ring(5)
        self.assertEqual(top.get_average_shortest_path_length(), 1.5)
        top = topology.ring(9)
        self.assertEqual(top.get_average_shortest_path_length(), 2.5)
        top = topology.custom()
        top.push_back()
        top.push_back()
        top.push_back()
        top.push_back()
        top.add_edge(0, 1)
        top.add_edge(1, 2)
        top.add_edge(2, 3)
        top.add_edge(3, 0)
        top.add_edge(0, 2)
        self.assertEqual(top.get_average_shortest_path_length(), 1.75)


def get_topology_test_suite():
    suite = _ut.TestSuite()
    suite.addTests(_ut.makeSuite(TopologyTests))
    return suite
