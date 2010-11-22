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

# This class will test that the C++ serialization methods (cpp_loads and cpp_dumps) are consistent.
class _serialization_test(_ut.TestCase):
	def __test_impl(self,types):
		for t in types:
			tmp1 = t()
			dump1 = tmp1.cpp_dumps()
			tmp2 = t()
			tmp2.cpp_loads(dump1)
			dump2 = tmp2.cpp_dumps()
			self.assertEqual(dump1,dump2)
	def test_problems(self):
		from PyGMO import problem
		types = filter(lambda t: not isinstance(t(),problem.base),[problem.__dict__[n] for n in filter(lambda n: not n.startswith('_') and not n == 'base',dir(problem))])
		self.__test_impl(types)
	def test_algorithms(self):
		from PyGMO import algorithm
		types = filter(lambda t: not isinstance(t(),algorithm.base),[algorithm.__dict__[n] for n in filter(lambda n: not n.startswith('_') and not n == 'base',dir(algorithm))])
		self.__test_impl(types)

# This class will stress the py_island class with highly concurrent simple evolutions.
class _py_island_torture_test(_ut.TestCase):
	def test_cpp(self):
		# Test with algo and prob implemented in C++.
		from PyGMO import py_island, archipelago, topology, algorithm, problem
		prob = problem.dejong(1)
		algo = algorithm.de(5)
		a = archipelago(topology.ring())
		for i in range(0,100):
			a.push_back(py_island(prob,algo,n = 20))
		a.evolve(10)
		a.join()
	def test_python(self):
		# Test with problem implemented in Python.
		from PyGMO import py_island, archipelago, topology, algorithm, problem
		prob = problem.py_test()
		algo = algorithm.de(5)
		a = archipelago(topology.ring())
		for i in range(0,100):
			a.push_back(py_island(prob,algo,n = 20))
		a.evolve(10)
		a.join()

def run_full_test_suite():
	"""
	Run the complete test suite for PyGMO.
	"""
	from PyGMO import test
	suite = _ut.TestLoader().loadTestsFromModule(test)
	_ut.TextTestRunner(verbosity=2).run(suite)
