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

def run_full_test_suite():
    """Run the complete test suite for PyGMO."""
    import sys
    from PyGMO import test
    from PyGMO.test._hypervolume_tests import get_hv_suite
    from PyGMO.test._topology_tests import get_topology_test_suite
    from PyGMO.test._archipelago_tests import get_archipelago_test_suite
    from PyGMO.test._serialization_tests import get_serialization_test_suite
    from PyGMO.test._island_torture_tests import get_island_torture_test_suite
    suite = _ut.TestLoader().loadTestsFromModule(test)

    # Add external suites explicitly
    suite.addTests(get_island_torture_test_suite())
    suite.addTests(get_serialization_test_suite())
    suite.addTests(get_hv_suite())
    suite.addTests(get_topology_test_suite())
    suite.addTests(get_archipelago_test_suite())

    successful = _ut.TextTestRunner(verbosity=2).run(suite).wasSuccessful()
    sys.exit(0 if successful else 1)

def run_island_torture_test_suite():
    """Run the island torture test suite."""
    from PyGMO.test._island_torture_tests import get_island_torture_test_suite
    _ut.TextTestRunner(verbosity=2).run(get_island_torture_test_suite())


def run_serialization_test_suite():
    """Run the serialization test suite."""
    from PyGMO.test._serialization_tests import get_serialization_test_suite
    _ut.TextTestRunner(verbosity=2).run(get_serialization_test_suite())


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
