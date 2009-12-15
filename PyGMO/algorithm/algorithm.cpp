/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

// 13/02/2008: Initial version by Francesco Biscani.

#include <boost/python/class.hpp>
#include <boost/python/module.hpp>

#include "../../src/GOclasses/algorithms/asa.h"
#include "../../src/GOclasses/algorithms/de.h"
#include "../../src/GOclasses/algorithms/mpso.h"
#include "../../src/GOclasses/algorithms/pso.h"
#include "../../src/GOclasses/algorithms/sga.h"
#include "../../src/GOclasses/algorithms/cs.h"
#include "../../src/GOclasses/algorithms/base.h"
#include "../../src/GOclasses/algorithms/ihs.h"
#include "../../src/GOclasses/algorithms/nm.h"
#include "../exceptions.h"
#include "../utils.h"

using namespace boost::python;
using namespace pagmo;

template <class Algorithm>
static inline class_<Algorithm,bases<algorithm::base> > algorithm_wrapper(const char *name, const char *descr)
{
	class_<Algorithm,bases<algorithm::base> > retval(name,descr,init<const Algorithm &>());
	retval.def("__copy__", &Py_copy_from_ctor<Algorithm>);
	retval.def("__repr__", &Py_repr_from_stream<Algorithm>);
	retval.add_property("id_name", &Algorithm::id_name, "Identification name.");
	retval.add_property("id_object", &Algorithm::id_object, "Object identification name.");
	return retval;
}

BOOST_PYTHON_MODULE(_algorithm) {
	// Translate exceptions for this module.
	translate_exceptions();

	// Expose base algorithm class.
	class_<algorithm::base, boost::noncopyable>("__algorithm::base", "Base GO algorithm", no_init);

	// Expose algorithms.
	algorithm_wrapper<algorithm::asa>("asa", "Simulated Annealing with adaptive neighbourhood algorithm.").def(init<int, const double &, const double &>()).def(init<int, const double&, const double&, const int, const int, const double>());
	algorithm_wrapper<algorithm::cs>("cs", "Compass search algorithm.").def(init<const double &, const double &, const double &>())
		.def(init<const double&>());
	algorithm_wrapper<algorithm::de>("de", "Differential evolution algorithm.").def(init<int, const double &, const double &, int>());
	algorithm_wrapper<algorithm::ihs>("ihs", "Improved harmony search algorithm.")
		.def(init<int, const double &, const double &, const double &, const double &, const double &>())
		.def(init<int>());
	algorithm_wrapper<algorithm::nm>("nm", "Nelder-Mead algorithm.")
		.def(init<int, const double &, const double &, const double &, const double &>())
		.def(init<int>());
	algorithm_wrapper<algorithm::mpso>("mpso", "MPSO algorithm.")
		.def(init<int, const double &, const double &, const double &, const double &, int>());
	algorithm_wrapper<algorithm::pso>("pso", "Particle swarm optimization algorithm.")
		.def(init<int, const double &, const double &, const double &, const double &>());
	algorithm_wrapper<algorithm::sga>("sga", "Simple genetic algorithm.")
		.def(init<int, const double &, const double &, int>())
		.def(init<int, const double &, const double &, int, double, int, int>());
}
