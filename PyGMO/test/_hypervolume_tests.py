from PyGMO import *
import unittest

class HVCtorTest(unittest.TestCase):

	def setUp(self):
		self.good_ps_2d = [[5,2],[3,4],[1,5]]

	def test_dimension_size(self):
		# fixed points dimensions should be reasonable
		self.assertRaises(ValueError, hypervolume, [[1,], [2,]])  # points of f_dim = 1
		self.assertRaises(ValueError, hypervolume, [[],[],])  # points of f_dim = 0
		self.assertRaises(ValueError, hypervolume, [[],])  # empty set of points
		self.assertRaises(ValueError, hypervolume, [[1,2,3], [3,3]])  # point dimensions are not equal
		self.assertRaises(TypeError, hypervolume, [])  # empty set of points #2
	
	def test_hypervolume_ctor_points(self):
		# test various possibilities of construction
		self.assertRaises(TypeError, hypervolume, "A string")  # something other than a list or a string
		self.assertRaises(TypeError, hypervolume, [[1,2,3], [2,3,4], [2,3,5]], "extra arg")  # good point definition with extra arg (hypervolume(ps, "something extra")
		self.assertRaises(TypeError, hypervolume, [[1,2,3], [2,3,4], "foo"])  # bad point
		self.assertRaises(TypeError, hypervolume, [[1,2,3], [2,3,4], [2,3,"bar"]])  # bad point value
		self.assertRaises(TypeError, hypervolume)  # skipping argument: hypervolume() should raise TypeError

		self.assertRaises(TypeError, hypervolume, foo=self.good_ps_2d)  # bad kwarg
		self.assertRaises(TypeError, hypervolume, self.good_ps_2d, foo="bar") # extra kwarg


class HVComputeTest(unittest.TestCase):

	def setUp(self):
		self.hv2d = hypervolume([[3,1],[2,2],[1,3]])

	def test_correct_out(self): 
		# simple 3D test
		hv = hypervolume([[1,1,1],[2,2,2,]])
		self.assertEqual(hv.compute(r = [3,3,3]), 8)

		# simple 2D test
		hv = hypervolume([[1,2],[2,1]])
		self.assertEqual(hv.compute(r = [3,3]), 3)

		# point on the border of refpoint (2D)
		hv = hypervolume([[1,2],[2,1]])
		self.assertEqual(hv.compute([2,2]), 0)

		# points on the border of refpoint (3D)
		hv = hypervolume([[1,2,1],[2,1,1]])
		self.assertEqual(hv.compute([2,2,2]), 0)

	
	def test_tuple_ctor(self):
		# test that hypervolume can be computed using a tuple as well
		hv = hypervolume(((1,1,1),(2,2,2,)))
		self.assertEqual(hv.compute(r = (3,3,3)), 8)

	def test_casting_float(self):
		# casting to float
		self.assertEqual(self.hv2d.compute(["4.0","4"]), 6)

	def test_refpoint_not_dom(self):
		# refpoint must be at least weakly dominated by every point (assuming minimization problem)
		hv = hypervolume([[1,3],[2,2], [3,1]])
		self.assertRaises(ValueError, hv.compute, [3,1])  # equal to some other point
		self.assertRaises(ValueError, hv.compute, [1.5,1.5])  # refpoint dominating some points
		self.assertRaises(ValueError, hv.compute, [0,0])  # refpoint dominating all points

	def test_kwargs(self):
		self.assertEqual(self.hv2d.compute(r=[3.5, 3.5]), 3.25)  # using kwarg 'r' correctly
		self.assertEqual(self.hv2d.compute(r=[3.5, 3.5], algorithm=hv_algorithm.native2d()), 3.25)  # using kwargs 'r', 'algorithm' correctly
	
		self.assertRaises(TypeError, self.hv2d.compute, refpoint=[4, 4])  # bad kwarg for reference point
		self.assertRaises(TypeError, self.hv2d.compute, r=[4, 4], foo="Something extra")  # we do not accept random kwargs
		self.assertRaises(TypeError, self.hv2d.compute, [4, 4], foo="Something extra")  # we do not accept random kwargs (as above but with ref point as arg)
		self.assertRaises(TypeError, self.hv2d.compute, [4, 4], hv_algorithm.native2d(), foo="Something extra")  # we do not accept random kwargs

	def test_kwargs_hv_algo(self):
		self.assertEqual(self.hv2d.compute(r=[4,4], algorithm=hv_algorithm.native2d()), 6)  # using kwargs correctly
		self.assertEqual(self.hv2d.compute(algorithm=hv_algorithm.native2d(), r=[4,4]), 6)  # using kwargs in reversed order
		self.assertEqual(self.hv2d.compute([4,4], algorithm=hv_algorithm.native2d()), 6)  # arg + kwarg
		self.assertEqual(self.hv2d.compute([4,4], hv_algorithm.native2d()), 6)  # arg + arg
		self.assertRaises(TypeError, self.hv2d.compute, algorithm=hv_algorithm.native2d())  # should use nadir point

	def test_bad_algo(self):
		self.assertRaises(ValueError, self.hv2d.compute, [4, 4], hv_algorithm.beume3d())  # 3d method to 2d problem
		self.assertRaises(ValueError, self.hv2d.compute, [4, 4], hv_algorithm.lebmeasure())  # lebmeasure to 2d problem

class HVLeastContribTest(unittest.TestCase):

	def setUp(self):
		self.r = [4,4]
		self.hv2d_eq_0 = hypervolume([[3,1],[2,2],[1,3]])  # all are equal (take first -> idx = 0)
		self.hv2d_eq_1 = hypervolume([[2.5,1],[2,2],[1,3]])  # last two are equal (take first -> idx = 1)
		self.hv2d_0 = hypervolume([[3.5,1],[2,2],[1,3]])
		self.hv2d_1 = hypervolume([[3,1],[2.5,2.5],[1,3]])
		self.hv2d_2 = hypervolume([[3,1],[2,2],[1,3.5]])
	
	def test_correct_out(self):
		self.assertEqual(self.hv2d_eq_0.least_contributor(r=self.r), 0)
		self.assertEqual(self.hv2d_eq_1.least_contributor(r=self.r), 1)
		self.assertEqual(self.hv2d_0.least_contributor(r=self.r), 0)
		self.assertEqual(self.hv2d_1.least_contributor(r=self.r), 1)
		self.assertEqual(self.hv2d_2.least_contributor(r=self.r), 2)
	
	def test_kwargs(self):
		self.assertEqual(self.hv2d_1.least_contributor(r=self.r), 1)  # using kwarg 'r' correctly
		self.assertEqual(self.hv2d_1.least_contributor(r=self.r, algorithm=hv_algorithm.native2d()), 1)  # using kwarg 'r' and 'algorithm' correctly

		self.assertRaises(TypeError, self.hv2d_0.least_contributor, refpoint=[4, 4])  # bad kwarg for reference point
		self.assertRaises(TypeError, self.hv2d_0.least_contributor, r=[4, 4], foo="Something extra")  # we do not accept random kwargs
		self.assertRaises(TypeError, self.hv2d_0.least_contributor, [4, 4], foo="Something extra")  # we do not accept random kwargs (as above but with ref point as arg)
		self.assertRaises(TypeError, self.hv2d_0.least_contributor, [4, 4], hv_algorithm.native2d(), foo="Something extra")  # we do not accept random kwargs

	def test_bad_algo(self):
		self.assertRaises(ValueError, self.hv2d_0.least_contributor, [4, 4], hv_algorithm.beume3d())  # 3d method to 2d problem
		self.assertRaises(ValueError, self.hv2d_0.least_contributor, [4, 4], hv_algorithm.lebmeasure())  # lebmeasure to 2d problem

class HVExclusiveTest(unittest.TestCase):
	def setUp(self):
		self.r = [4,4]
		self.hv2d = hypervolume([[3,1],[2,2],[1,3]])  # all are equal (take first -> idx = 0)
		self.hv2d_2 = hypervolume([[3.1,1],[2,2],[1,3]])  # all are equal (take first -> idx = 0)
	
	def test_correct_out(self):
		self.assertEqual(self.hv2d.exclusive(p_idx=0, r=self.r), 1)
		self.assertEqual(self.hv2d.exclusive(p_idx=1, r=self.r), 1)
		self.assertEqual(self.hv2d.exclusive(p_idx=2, r=self.r), 1)
		self.assertTrue(abs(self.hv2d_2.exclusive(p_idx=0, r=self.r) - 0.9) < 0.00000001)
	
	def test_kwargs(self):
		self.assertRaises(TypeError,self.hv2d.exclusive, 0)  # no refpoint
		self.assertRaises(TypeError,self.hv2d.exclusive, p_idx=0)  # no refpoint
		self.assertRaises(TypeError,self.hv2d.exclusive, p_idx=0, hv_algorithm=hv_algorithm.wfg())  # no refpoint
		self.assertEqual(self.hv2d.exclusive(0, self.r), 1)  # using arg
		self.assertEqual(self.hv2d.exclusive(0, r=self.r), 1)  # using kwarg 'r' correctly
		self.assertEqual(self.hv2d.exclusive(p_idx=0, r=self.r), 1)  # using kwarg 'r' correctly

		self.assertRaises(TypeError, self.hv2d.exclusive, p_idx=0, algorithm=hv_algorithm.native2d())  # no refpoint
		self.assertEqual(self.hv2d.exclusive(0, self.r, hv_algorithm.native2d()), 1)  # all args
		self.assertEqual(self.hv2d.exclusive(0, self.r, algorithm=hv_algorithm.native2d()), 1)  # last kwarg
		self.assertEqual(self.hv2d.exclusive(p_idx=0, r=self.r, algorithm=hv_algorithm.native2d()), 1)  # all kwargs
		self.assertEqual(self.hv2d.exclusive(algorithm=hv_algorithm.native2d(), r=self.r, p_idx=0), 1)  # all kwargs in reverse

		self.assertRaises(TypeError, self.hv2d.exclusive, 0, refpoint=[4, 4])  # bad kwarg for reference point
		self.assertRaises(TypeError, self.hv2d.exclusive, 0, r=[4, 4], foo="Something extra")  # we do not accept random kwargs
		self.assertRaises(TypeError, self.hv2d.exclusive, 0, [4, 4], foo="Something extra")  # we do not accept random kwargs (as above but with ref point as arg)
		self.assertRaises(TypeError, self.hv2d.exclusive, 0, [4, 4], hv_algorithm.native2d(), foo="Something extra")  # we do not accept random kwargs
		self.assertRaises(TypeError, self.hv2d.exclusive, p_idx=0, r=self.r, algorithm=hv_algorithm.native2d(), foo="Something extra")
		self.assertRaises(TypeError, self.hv2d.exclusive, r=self.r, algorithm=hv_algorithm.native2d()) # no p_idx 
	
	def test_p_idx(self):
		self.assertRaises(TypeError, self.hv2d.exclusive, -1, self.r)  # negative idx
		self.assertRaises(ValueError, self.hv2d.exclusive, 100, self.r)  # large
		self.assertRaises(TypeError, self.hv2d.exclusive, "not an int", self.r)  # large
		self.assertRaises(TypeError, self.hv2d.exclusive)  # p_idx not provided

	def test_bad_algo(self):
		self.assertRaises(ValueError, self.hv2d.exclusive, 0, [4, 4], hv_algorithm.beume3d())  # 3d method to 2d problem
		self.assertRaises(ValueError, self.hv2d.exclusive, 0, [4, 4], hv_algorithm.lebmeasure())  # lebmeasure to 2d problem

class HVNadirPointTest(unittest.TestCase):

	def setUp(self):
		self.hv2d = hypervolume([[3,1],[2,2],[1,3]])

	def test_nadir_point(self):
		self.assertEqual(tuple(self.hv2d.get_nadir_point()), (4,4))  # default nadir point
		self.assertEqual(tuple(self.hv2d.get_nadir_point(5.0)), (8,8))  # custom nadir point
		self.assertEqual(tuple(self.hv2d.get_nadir_point(0.0)), (3,3))  # nadir point with eps=0.0
		self.assertEqual(tuple(self.hv2d.get_nadir_point(-0.0)), (3,3))  # nadir point with eps=-0.0 is ok
		self.assertEqual(tuple(self.hv2d.get_nadir_point(eps=5.0)), (8,8))  # custom nadir point with 'eps' kwarg
		self.assertRaises(ValueError, self.hv2d.get_nadir_point, -0.0000001)  # nadir point with negative eps
		self.assertRaises(TypeError, self.hv2d.get_nadir_point, eps="foo")  # bad kwarg
		self.assertRaises(TypeError, self.hv2d.get_nadir_point, "foo") # bad arg
		self.assertRaises(TypeError, self.hv2d.get_nadir_point, epsilon=1.0)  # bad kwarg name

def get_hv_suite():
	suite = unittest.TestSuite()
	suite.addTests(unittest.makeSuite(HVCtorTest))
	suite.addTests(unittest.makeSuite(HVComputeTest))
	suite.addTests(unittest.makeSuite(HVLeastContribTest))
	suite.addTests(unittest.makeSuite(HVExclusiveTest))
	suite.addTests(unittest.makeSuite(HVNadirPointTest))
	return suite
