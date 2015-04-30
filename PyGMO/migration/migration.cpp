/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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

// Workaround for http://mail.python.org/pipermail/new-bugs-announce/2011-March/010395.html
#ifdef _WIN32
#include <cmath>
#endif
 
#include <boost/python/class.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/module.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/utility.hpp>

#include "../../src/migration.h"
#include "../utils.h"

using namespace boost::python;
using namespace pagmo;

template <class MSPolicy>
static inline class_<MSPolicy,bases<migration::base_s_policy> > migration_s_policy_wrapper(const char *name, const char *descr)
{
	class_<MSPolicy,bases<migration::base_s_policy> > retval(name,descr,init<const MSPolicy &>());
	retval.def(init<>());
	retval.def("__copy__", &Py_copy_from_ctor<MSPolicy>);
	retval.def("__deepcopy__", &Py_deepcopy_from_ctor<MSPolicy>);
	retval.def("__repr__", &migration::base::human_readable);
	retval.def_pickle(python_class_pickle_suite<MSPolicy>());
	retval.def("cpp_loads", &py_cpp_loads<MSPolicy>);
	retval.def("cpp_dumps", &py_cpp_dumps<MSPolicy>);
	return retval;
}

template <class MRPolicy>
static inline class_<MRPolicy,bases<migration::base_r_policy> > migration_r_policy_wrapper(const char *name, const char *descr)
{
	class_<MRPolicy,bases<migration::base_r_policy> > retval(name,descr,init<const MRPolicy &>());
	retval.def(init<>());
	retval.def("__copy__", &Py_copy_from_ctor<MRPolicy>);
	retval.def("__deepcopy__", &Py_deepcopy_from_ctor<MRPolicy>);
	retval.def("__repr__", &migration::base::human_readable);
	retval.def_pickle(python_class_pickle_suite<MRPolicy>());
	retval.def("cpp_loads", &py_cpp_loads<MRPolicy>);
	retval.def("cpp_dumps", &py_cpp_dumps<MRPolicy>);
	return retval;
}

BOOST_PYTHON_MODULE(_migration) {
	common_module_init();

	// Migration rate type enum.
	enum_<migration::rate_type>("rate_type")
		.value("absolute",migration::absolute)
		.value("fractional",migration::fractional);

	// Expose migration selection policies.

	// Base.
	class_<migration::base_s_policy,boost::noncopyable>("_base_s_policy",no_init);

	// Best selection policy.
	migration_s_policy_wrapper<migration::best_s_policy>("best_s_policy","Best migration selection policy.")
		.def(init<optional<const double &, migration::rate_type> >());

	// Random selection policy.
	migration_s_policy_wrapper<migration::random_s_policy>("random_s_policy","Selection policy for random individuals.")
		.def(init<optional<const double &, migration::rate_type> >());

	// Best Kill selection policy.
	migration_s_policy_wrapper<migration::best_kill_s_policy>("best_kill_s_policy","Best Kill migration selection policy.")
		.def(init<optional<const double &, migration::rate_type> >());

	// Hypervolume greedy selection policy
	migration_s_policy_wrapper<migration::hv_greedy_s_policy>("hv_greedy_s_policy","Hypervolume Greedy migration selection policy.")
		.def(init<optional<const double &, migration::rate_type, const double> >());

	// Hypervolume best selection policy
	migration_s_policy_wrapper<migration::hv_best_s_policy>("hv_best_s_policy","Hypervolume Best migration selection policy.")
		.def(init<optional<const double &, migration::rate_type, const double> >());

	// Expose migration replacement policies.	

	// Base.
	class_<migration::base_r_policy,boost::noncopyable>("_base_r_policy",no_init);

	// Fair replacement policy.
	migration_r_policy_wrapper<migration::fair_r_policy>("fair_r_policy","Fair migration replacement policy.")
		.def(init<optional<const double &, migration::rate_type> >());

	// Random replacement policy.
	migration_r_policy_wrapper<migration::random_r_policy>("random_r_policy","Random migration replacement policy.")
		.def(init<optional<const double &, migration::rate_type> >());

	// Worst replacement policy.
	migration_r_policy_wrapper<migration::worst_r_policy>("worst_r_policy","Worst migration replacement policy.")
		.def(init<optional<const double &, migration::rate_type> >());

	// Hypervolume greedy replacement policy.
	migration_r_policy_wrapper<migration::hv_greedy_r_policy>("hv_greedy_r_policy","Hypervolume greedy migration replacement policy.")
		.def(init<optional<const double &, migration::rate_type, const double> >());

	// Hypervolume fair replacement policy.
	migration_r_policy_wrapper<migration::hv_fair_r_policy>("hv_fair_r_policy","Hypervolume fair migration replacement policy.")
		.def(init<optional<const double &, migration::rate_type, const double> >());

	// Register to_python conversion from smart pointer.
	register_ptr_to_python<migration::base_s_policy_ptr>();
	register_ptr_to_python<migration::base_r_policy_ptr>();
}
