from PyGMO import population, problem
from PyGMO.util import hypervolume, hv_algorithm
import unittest


class HVCtorTest(unittest.TestCase):

    def setUp(self):
        self.good_ps_2d = [[5, 2], [3, 4], [1, 5]]

    def test_dimension_size(self):
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


class HVFlagsTest(unittest.TestCase):

    def test_gettersetter(self):
        hv = hypervolume([[1, 2, 3], [4, 5, 6]])
        self.assertTrue(hv.get_verify())
        self.assertTrue(hv.get_copy_points())
        hv.set_verify(False)
        hv.set_copy_points(False)
        self.assertTrue(hv.get_verify() is False)
        self.assertTrue(hv.get_copy_points() is False)


class HVComputeTest(unittest.TestCase):

    def setUp(self):
        self.hv2d = hypervolume([[3, 1], [2, 2], [1, 3]])

    def test_correct_out(self):
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
        hv = hypervolume([[1, 1, 1, 1], [2, 2, 2, 2]])
        self.assertEqual(
            hv.compute(r=(3, 3, 3, 3), algorithm=hv_algorithm.hv4d()), 16)
        self.assertEqual(
            hv.compute(r=(3, 3, 3, 3), algorithm=hv_algorithm.fpl()), 16)

    def test4d_edge(self):
        hv = hypervolume([[1, 1, 1, 3], [2, 2, 2, 3]])
        self.assertEqual(
            hv.compute(r=(3, 3, 3, 3), algorithm=hv_algorithm.hv4d()), 0)
        self.assertEqual(
            hv.compute(r=(3, 3, 3, 3), algorithm=hv_algorithm.fpl()), 0)

    def test4d_duplicate(self):
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
        hv = hypervolume(((1, 1, 1), (2, 2, 2,)))
        self.assertEqual(hv.compute(r=(3, 3, 3)), 8)

    def test_casting_float(self):
        # casting to float
        self.assertEqual(self.hv2d.compute(["4.0", "4"]), 6)

    def test_refpoint_not_dom(self):
        # refpoint must be at least weakly dominated by every point (assuming
        # minimization problem)
        hv = hypervolume([[1, 3], [2, 2], [3, 1]])
        # equal to some other point
        self.assertRaises(ValueError, hv.compute, [3, 1])
        # refpoint dominating some points
        self.assertRaises(ValueError, hv.compute, [1.5, 1.5])
        # refpoint dominating all points
        self.assertRaises(ValueError, hv.compute, [0, 0])

    def test_kwargs(self):
        # using kwarg 'r' correctly
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
        # 3d method to 2d problem
        self.assertRaises(
            ValueError, self.hv2d.compute, [4, 4], hv_algorithm.hv3d())


class HVContributionsTest(unittest.TestCase):

    def assertContribs(self, S, R, ans):
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


class HVLeastContribTest(unittest.TestCase):

    def setUp(self):
        self.r = [4, 4]
        self.hv2d_eq_0 = hypervolume([[3, 1], [2, 2], [1, 3]])  # LC in [0,1,2]
        self.hv2d_eq_1 = hypervolume([[2.5, 1], [2, 2], [1, 3]])  # LC = 1
        self.hv2d_0 = hypervolume([[3.5, 1], [2, 2], [1, 3]])
        self.hv2d_1 = hypervolume([[3, 1], [2.5, 2.5], [1, 3]])
        self.hv2d_2 = hypervolume([[3, 1], [2, 2], [1, 3.5]])

    def test_correct_out(self):
        self.assertTrue(
            self.hv2d_eq_0.least_contributor(r=self.r) in [0, 1, 2])

        self.assertEqual(self.hv2d_eq_1.least_contributor(r=self.r), 1)
        self.assertEqual(self.hv2d_0.least_contributor(r=self.r), 0)
        self.assertEqual(self.hv2d_1.least_contributor(r=self.r), 1)
        self.assertEqual(self.hv2d_2.least_contributor(r=self.r), 2)

    def test_kwargs(self):
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
        self.assertRaises(
            ValueError, self.hv2d_0.least_contributor, [
                4, 4], hv_algorithm.hv3d())  # 3d method to 2d problem


class HVExclusiveTest(unittest.TestCase):

    def setUp(self):
        self.r = [4, 4]
        # all are equal (take first -> idx = 0)
        self.hv2d = hypervolume([[3, 1], [2, 2], [1, 3]])
        # all are equal (take first -> idx = 0)
        self.hv2d_2 = hypervolume([[3.1, 1], [2, 2], [1, 3]])

    def test_correct_out(self):
        self.assertEqual(self.hv2d.exclusive(p_idx=0, r=self.r), 1)
        self.assertEqual(self.hv2d.exclusive(p_idx=1, r=self.r), 1)
        self.assertEqual(self.hv2d.exclusive(p_idx=2, r=self.r), 1)
        self.assertTrue(
            abs(self.hv2d_2.exclusive(p_idx=0, r=self.r) - 0.9) < 0.00000001)

    def test_kwargs(self):
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
        # negative idx
        self.assertRaises(TypeError, self.hv2d.exclusive, -1, self.r)
        self.assertRaises(
            ValueError, self.hv2d.exclusive, 100, self.r)  # large
        self.assertRaises(
            TypeError, self.hv2d.exclusive, "not an int", self.r)  # large
        self.assertRaises(TypeError, self.hv2d.exclusive)  # p_idx not provided

    def test_bad_algo(self):
        # 3d method to 2d problem
        self.assertRaises(
            ValueError, self.hv2d.exclusive, 0, [4, 4], hv_algorithm.hv3d())


class HVNadirPointTest(unittest.TestCase):

    def setUp(self):
        self.hv2d = hypervolume([[3, 1], [2, 2], [1, 3]])

    def test_nadir_point(self):
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


class HVAlgorithms(unittest.TestCase):

    def setUp(self):
        self.ps4d_0 = [[1., ] * 4, ]
        self.ps3d_0 = [[1., ] * 3, ]
        self.r4d_0 = [2, ] * 4
        self.r3d_0 = [2, ] * 3

    def test_wfg(self):
        # stop_dimension = 0
        self.assertRaises(ValueError, hv_algorithm.wfg, stop_dimension=0)
        # stop_dimension = 1
        self.assertRaises(ValueError, hv_algorithm.wfg, stop_dimension=1)

    def test_hv4d(self):
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
    suite = unittest.TestSuite()
    suite.addTests(unittest.makeSuite(HVCtorTest))
    suite.addTests(unittest.makeSuite(HVComputeTest))
    suite.addTests(unittest.makeSuite(HVLeastContribTest))
    suite.addTests(unittest.makeSuite(HVExclusiveTest))
    suite.addTests(unittest.makeSuite(HVNadirPointTest))
    suite.addTests(unittest.makeSuite(HVAlgorithms))
    suite.addTests(unittest.makeSuite(HVFlagsTest))
    suite.addTests(unittest.makeSuite(HVContributionsTest))
    return suite
