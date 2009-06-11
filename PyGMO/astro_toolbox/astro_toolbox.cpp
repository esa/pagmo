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

// 05/03/2009: Initial version by Francesco Biscani.

#include <boost/python.hpp>
#include <vector>

#include "../../src/AstroToolbox/Astro_Functions.h"
#include "../../src/AstroToolbox/Lambert.h"
#include "../../src/AstroToolbox/propagateKEP.h"
#include "../../src/keplerian_toolbox/kstate.h"
#include "../boost_python_container_conversions.h"
#include "../exceptions.h"

using namespace boost::python;
using namespace keplerian_toolbox;

struct lambert_result {
	lambert_result():v0(size_t(3)),v1(size_t(3)),a(0),p(0),theta(0),it(0) {}
	std::vector<double> v0;
	std::vector<double> v1;
	double				a;
	double				p;
	double				theta;
	int 				it;
};

static inline void Py_lambertI(const std::vector<double> &r0, const std::vector<double> &r1, const double &T, const double &mu, bool lw,
	lambert_result &lr)
{
	if (r0.size() != 3 || r1.size() != 3 || lr.v0.size() != 3 || lr.v1.size() != 3) {
		P_EX_THROW(value_error,"the size of all input/output position/velocity vectors must be 3");
	}
	LambertI(&r0[0],&r1[0],T,mu,lw,&lr.v0[0],&lr.v1[0],lr.a,lr.p,lr.theta,lr.it);
}

static inline void Py_propagate_kep(const std::vector<double> &r0, const std::vector<double> &v0, const double &t, const double &mu,
	std::vector<double> &r1, std::vector<double> &v1)
{
	if (r0.size() != 3 || r1.size() != 3 || v0.size() != 3 || v1.size() != 3) {
		P_EX_THROW(value_error,"the size of all input/output position/velocity vectors must be 3");
	}
	propagateKEP(&r0[0],&v0[0],t,mu,&r1[0],&v1[0]);
}

BOOST_PYTHON_MODULE(_astro_toolbox) {
	to_tuple_mapping<array_d3>();
	from_python_sequence<std::vector<double>,variable_capacity_policy>();

	// Translate exceptions for this module.
	translate_exceptions();

	class_<lambert_result>("__lambert_result",init<>())
		.def_readonly("v0", &lambert_result::v0)
		.def_readonly("v1", &lambert_result::v1)
		.def_readonly("a", &lambert_result::a)
		.def_readonly("p", &lambert_result::p)
		.def_readonly("theta", &lambert_result::theta)
		.def_readonly("it", &lambert_result::it);
	def("__lambertI", &Py_lambertI);
	def("__propagate_kep", &Py_propagate_kep);
	def("M2E", &Mean2Eccentric, "Convert mean anomaly M into eccentric anomaly E by solving Kepler's equation.");

	class_<kstate> ks("kstate",init<>());
	ks.def(init<const std::vector<double> &, const std::vector<double> &, const double &, const double &>());
	ks.def(repr(self));
	ks.add_property("r",make_function(&kstate::get_r, return_value_policy<return_by_value>()), &kstate::set_r<std::vector<double> >);
	ks.add_property("v",make_function(&kstate::get_v, return_value_policy<return_by_value>()), &kstate::set_v<std::vector<double> >);
	ks.def("propagate", &kstate::propagate, return_value_policy<return_by_value>());
}
