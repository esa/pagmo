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
	def __test_impl(self,types):
		for t in types:
			tmp = t()
			dump1 = tmp.cpp_dumps()
			tmp.cpp_loads(dump1)
			dump2 = tmp.cpp_dumps()
			self.assertEqual(dump1,dump2)
	def test_problems(self):
		from PyGMO import problem
		types = filter(lambda t: not isinstance(t(),problem.base),[problem.__dict__[n] for n in filter(lambda n: not n.startswith('_') and not n == 'base',dir(problem))])
		self.__test_impl(types)
	def test_algorithms(self):
		from PyGMO import algorithm
		types = filter(lambda t: not isinstance(t(),algorithm.base),[algorithm.__dict__[n] for n in filter(lambda n: not n.startswith('_') and not n == 'base',dir(algorithm))])
		self.__test_impl(types)

def run_full_test_suite():
	suite = _ut.TestLoader().loadTestsFromTestCase(_serialization_test)
	_ut.TextTestRunner(verbosity=2).run(suite)
