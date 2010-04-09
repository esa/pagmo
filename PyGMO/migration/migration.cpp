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

#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/register_ptr_to_python.hpp>

#include "../../src/migration.h"
#include "../exceptions.h"

using namespace boost::python;
using namespace pagmo;

template <class MPolicy>
static inline class_<MPolicy> migration_policy_wrapper(const char *name, const char *descr)
{
	class_<MPolicy> retval(name,descr,init<const MPolicy &>());
	retval.def("__copy__", &MPolicy::clone);
	return retval;
}

BOOST_PYTHON_MODULE(_migration) {
	// Translate exceptions for this module.
	translate_exceptions();

	// Expose migration selection policies.

	// Best selection policy.
	migration_policy_wrapper<migration::best_s_policy>("best_s_policy","Best migration selection policy.")
		.def(init<>());

	// Expose migration replacement policies.	

	// Fair replacement policy.
	migration_policy_wrapper<migration::fair_r_policy>("fair_r_policy","Fair migration replacement policy.")
		.def(init<>());

	// Random replacement policy.
	migration_policy_wrapper<migration::random_r_policy>("random_r_policy","Random migration replacement policy.")
		.def(init<>());

	// Register to_python conversion from smart pointer.
	register_ptr_to_python<migration::base_s_policy_ptr>();
	register_ptr_to_python<migration::base_r_policy_ptr>();
}
