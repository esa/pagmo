/*****************************************************************************
 *   Copyright (C) 2008, 2009 Advanced Concepts Team (European Space Agency) *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
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
#include <boost_python_p_exceptions.h>

#include "../../src/GOclasses/algorithms/ASA.h"
#include "../../src/GOclasses/algorithms/DE.h"
#include "../../src/GOclasses/algorithms/MPSO.h"
#include "../../src/GOclasses/algorithms/PSO.h"
#include "../../src/GOclasses/algorithms/SGA.h"
#include "../../src/GOclasses/algorithms/CS.h"
#include "../../src/GOclasses/algorithms/go_algorithm.h"
#include "../../src/GOclasses/algorithms/ihs_algorithm.h"
#include "../../src/GOclasses/algorithms/nm_algorithm.h"
#include "../utils.h"

using namespace boost::python;

template <class Algorithm>
static inline class_<Algorithm,bases<go_algorithm> > algorithm_wrapper(const char *name, const char *descr)
{
	class_<Algorithm,bases<go_algorithm> > retval(name,descr,init<const Algorithm &>());
	retval.def("__copy__", &Py_copy_from_ctor<Algorithm>);
	retval.def("__repr__", &Py_repr_from_stream<Algorithm>);
	retval.add_property("id_name", &Algorithm::id_name, "Identification name.");
	retval.add_property("id_object", &Algorithm::id_object, "Object identification name.");
	return retval;
}

BOOST_PYTHON_MODULE(_algorithm) {
	// Translate exceptions for this module.
	translate_p_exceptions();

	// Expose base algorithm class.
	class_<go_algorithm, boost::noncopyable>("__go_algorithm", "Base GO algorithm", no_init);

	// Expose algorithms.
        algorithm_wrapper<ASAalgorithm>("asa", "Simulated Annealing with adaptive neighbourhood algorithm.").def(init<int, const double &, const double &>());
        algorithm_wrapper<CSalgorithm>("cs", "Compass search algorithm.").def(init<const double &, const double &, const double &>());
        algorithm_wrapper<DEalgorithm>("de", "Differential evolution algorithm.").def(init<int, const double &, const double &, int>());
	algorithm_wrapper<ihs_algorithm>("ihs", "IHS algorithm.")
		.def(init<int, const double &, const double &, const double &, const double &, const double &>());
	algorithm_wrapper<nm_algorithm>("nm", "Nelder-Mead algorithm.")
		.def(init<int, const double &, const double &, const double &, const double &>());
	algorithm_wrapper<MPSOalgorithm>("mpso", "MPSO algorithm.")
		.def(init<int, const double &, const double &, const double &, const double &, int>());
        algorithm_wrapper<PSOalgorithm>("pso", "Particle swarm optimization algorithm.")
		.def(init<int, const double &, const double &, const double &, const double &>());
	algorithm_wrapper<SGAalgorithm>("sga", "Simple genetic algorithm.")
		.def(init<int, const double &, const double &, int>());
}
