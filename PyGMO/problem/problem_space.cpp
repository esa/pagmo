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
 
#include <Python.h>
#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/module.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/utility.hpp>
#include <cstddef>
#include <string>

#include "../../src/exceptions.h"
#include "../../src/problems.h"
#include "../../src/serialization.h"
#include "../../src/types.h"
#include "../../src/config.h"
#include "../utils.h"
#include "python_base.h"
#include "python_base_stochastic.h"

#include <keplerian_toolbox/planet/base.h>
#include "../../src/problem/tsp_ds.h"


using namespace boost::python;
using namespace pagmo;



// wrapper for the find_subsequence method of tsp_ads
static inline tuple find_subsequence_wrapper_ds(const problem::tsp_ds& p, const decision_vector& tour, const bool static_computations)
{
	double retval_p, retval_l;
	decision_vector::size_type retval_it_l, retval_it_r;
    p.find_subsequence(tour,retval_p,retval_l,retval_it_l,retval_it_r, static_computations);
	return boost::python::make_tuple(retval_p,retval_l,retval_it_l,retval_it_r);
}


// Wrapper to expose problems.
template <class Problem>
static inline class_<Problem,bases<problem::base> > problem_wrapper(const char *name, const char *descr)
{
	class_<Problem,bases<problem::base> > retval(name,descr,init<const Problem &>());
	retval.def(init<>());
	retval.def("__copy__", &Py_copy_from_ctor<Problem>);
	retval.def("__deepcopy__", &Py_deepcopy_from_ctor<Problem>);
	retval.def_pickle(python_class_pickle_suite<Problem>());
	retval.def("cpp_loads", &py_cpp_loads<Problem>);
	retval.def("cpp_dumps", &py_cpp_dumps<Problem>);
	return retval;
}

// Wrapper to expose stochastic problems.
template <class Problem>
static inline class_<Problem,bases<problem::base>,bases<problem::base_stochastic> > stochastic_problem_wrapper(const char *name, const char *descr)
{
	class_<Problem,bases<problem::base>,bases<problem::base_stochastic> > retval(name,descr,init<const Problem &>());
	retval.def(init<>());
	retval.def("__copy__", &Py_copy_from_ctor<Problem>);
	retval.def("__deepcopy__", &Py_deepcopy_from_ctor<Problem>);
	retval.def_pickle(python_class_pickle_suite<Problem>());
	retval.def("cpp_loads", &py_cpp_loads<Problem>);
	retval.def("cpp_dumps", &py_cpp_dumps<Problem>);
	retval.add_property("seed",&problem::base_stochastic::get_seed,&problem::base_stochastic::set_seed,
		"Random seed used in the objective function evaluation.");
	return retval;
}


// Wrapper to expose TSP problems.
template <class Problem>
static inline class_<Problem,bases<problem::base>,bases<problem::base_tsp> > tsp_problem_wrapper(const char *name, const char *descr)
{
	class_<Problem,bases<problem::base>,bases<problem::base_tsp> > retval(name,descr,init<const Problem &>());
	retval.def(init<>());
	retval.def("__copy__", &Py_copy_from_ctor<Problem>);
	retval.def("__deepcopy__", &Py_deepcopy_from_ctor<Problem>);
	retval.def_pickle(python_class_pickle_suite<Problem>());
	retval.def("cpp_loads", &py_cpp_loads<Problem>);
	retval.def("cpp_dumps", &py_cpp_dumps<Problem>);
	retval.def("full2cities", &problem::base_tsp::full2cities);
	retval.def("cities2full", &problem::base_tsp::cities2full);
	retval.def("randomkeys2cities", &problem::base_tsp::randomkeys2cities);
	retval.def("cities2randomkeys", &problem::base_tsp::cities2randomkeys);
	retval.add_property("encoding", &problem::base_tsp::get_encoding);
	retval.add_property("n_cities", &problem::base_tsp::get_n_cities);
	return retval;
}


BOOST_PYTHON_MODULE(_problem_space) {
	common_module_init();


// Asteroid Sample Return (also used fot human missions to asteroids)
//	problem_wrapper<problem::sample_return>("sample_return","Asteroid sample return problem.")
//		.def(init<const ::kep_toolbox::planet &, optional<const double &> >())
//		.def("get_delta_v",&problem::sample_return::get_delta_v);

	// Earth-planet problem.
	//problem_wrapper<problem::earth_planet>("earth_planet","Earth-planet low-thrust problem.")
		//.def(init< optional<int, std::string, const double &> >());

	// GTOC5 launch.
	//problem_wrapper<problem::gtoc5_launch>("gtoc5_launch","GTOC5 launch phase.")
		//.def(init< optional<int, int, problem::gtoc5_launch::objective, const double &> >());

	//enum_<problem::gtoc5_launch::objective>("gtoc5_launch_objective")
		//.value("MASS",problem::gtoc5_launch::MASS)
		//.value("TIME",problem::gtoc5_launch::TIME);

	// GTOC5 randez-vouz.
	//problem_wrapper<problem::gtoc5_rendezvous>("gtoc5_rendezvous","GTOC5 rendezvous phase.")
		//.def(init< optional<int, int, int, const double &, const double &, const double &> >());

	//enum_<problem::gtoc5_flyby::objective>("gtoc5_flyby_objective")
		//.value("MASS",problem::gtoc5_flyby::MASS)
		//.value("FINAL_EPOCH",problem::gtoc5_flyby::FINAL_EPOCH)
		//.value("TIME",problem::gtoc5_flyby::TIME);

	// GTOC5 flyby.
	//problem_wrapper<problem::gtoc5_flyby>("gtoc5_flyby","GTOC5 flyby phase.")
		//.def(init< optional<int, int, int, int, const double &, const double &, problem::gtoc5_flyby::objective, const double &, const double &> >());

	// GTOC5 self flyby.
	//problem_wrapper<problem::gtoc5_self_flyby>("gtoc5_self_flyby","GTOC5 self flyby phase.")
		//.def(init< optional<int, int, const double &, const double &, const double &> >());

	// GTOC1 problem.
	problem_wrapper<problem::gtoc_1>("gtoc_1","GTOC 1 problem (chemical approximation).");

	// GTOC2 problem.
	problem_wrapper<problem::gtoc_2>("gtoc_2","GTOC 2 problem (LT model).")
		.def(init<optional<int,int,int,int,int,problem::gtoc_2::objective> >());

	// GTOC2's objectives enum.
	enum_<problem::gtoc_2::objective>("_gtoc_2_objective")
		.value("MASS",problem::gtoc_2::MASS)
		.value("TIME",problem::gtoc_2::TIME)
		.value("MASS_TIME",problem::gtoc_2::MASS_TIME);

	// Laplace problem.
	problem_wrapper<problem::laplace>("laplace","Laplace problem.")
		.def(init< const std::vector<int> &>());

	// Cassini 1.
	problem_wrapper<problem::cassini_1>("cassini_1","Cassini 1 interplanetary trajectory problem.")
		.def(init<optional<unsigned int> >());
	// Messenger full.
	problem_wrapper<problem::messenger_full>("messenger_full","Full Messenger problem.");

	// Cassini 2.
	problem_wrapper<problem::cassini_2>("cassini_2","Cassini 2 interplanetary trajectory problem.");

	// Rosetta problem.
	problem_wrapper<problem::rosetta>("rosetta","Rosetta problem.");

	// Sagas problem.
	problem_wrapper<problem::sagas>("sagas","Sagas problem.");

	// Tandem.
	problem_wrapper<problem::tandem>("tandem","Tandem problem.")
		.def(init< optional<int, double> >())
		.def("pretty", &problem::tandem::pretty);

	problem_wrapper<problem::mga_1dsm_alpha>("mga_1dsm_alpha", "A Multiple Gravity Assist with 1 Deep Space Manouvre problem")
		.def(init< optional<std::vector<kep_toolbox::planet::planet_ptr>, kep_toolbox::epoch, kep_toolbox::epoch, double, double, double, double, bool, bool, bool> >())
		.def("pretty", &problem::mga_1dsm_alpha::pretty)
		.def("set_tof", &problem::mga_1dsm_alpha::set_tof)
		.def("get_tof", &problem::mga_1dsm_alpha::get_tof)
		.def("set_launch_window", &problem::mga_1dsm_alpha::set_launch_window)
		.def("set_vinf", &problem::mga_1dsm_alpha::set_vinf)
		.def("get_sequence", &problem::mga_1dsm_alpha::get_sequence);

	problem_wrapper<problem::mga_1dsm_tof>("mga_1dsm_tof", "A Multiple Gravity Assist with 1 Deep Space Manouvre problem")
		.def(init< optional<std::vector<kep_toolbox::planet::planet_ptr>, kep_toolbox::epoch, kep_toolbox::epoch, std::vector<boost::array<double,2> >, double, double, bool, bool, bool> >())
		.def("pretty", &problem::mga_1dsm_tof::pretty, (arg("x"),arg("extended_output") = false))
		.def("set_tof", &problem::mga_1dsm_tof::set_tof)
		.def("get_tof", &problem::mga_1dsm_tof::get_tof)
		.def("set_launch_window", &problem::mga_1dsm_tof::set_launch_window)
		.def("set_vinf", &problem::mga_1dsm_tof::set_vinf)
		.def("get_sequence", &problem::mga_1dsm_tof::get_sequence);
		
	problem_wrapper<problem::mga_incipit>("mga_incipit", "Jupiter capture problem from the first part of gtoc6")
		.def(init< optional< std::vector<kep_toolbox::planet::planet_ptr>, kep_toolbox::epoch, kep_toolbox::epoch, std::vector<std::vector<double> > > >())
		.def("pretty", &problem::mga_incipit::pretty)
		.def("get_sequence", &problem::mga_incipit::get_sequence)
		.add_property("tof",make_function(&problem::mga_incipit::get_tof, return_value_policy<copy_const_reference>()), &problem::mga_incipit::set_tof,"bound on the times of flight for the different legs");

	problem_wrapper<problem::mga_incipit_cstrs>("mga_incipit_cstrs", "Jupiter capture problem from the first part of gtoc6 (constrained version)")
		.def(init< optional< std::vector<kep_toolbox::planet::planet_ptr>, kep_toolbox::epoch, kep_toolbox::epoch, std::vector<std::vector<double> >, const double, const double > >())
		.def("pretty", &problem::mga_incipit_cstrs::pretty)
		.def("get_sequence", &problem::mga_incipit_cstrs::get_sequence)
		.add_property("tof",make_function(&problem::mga_incipit_cstrs::get_tof, return_value_policy<copy_const_reference>()), &problem::mga_incipit_cstrs::set_tof,"bound on the times of flight for the different legs");

	problem_wrapper<problem::mga_part>("mga_part", "A part of the Jupiter moon tour from gtoc6")
		.def(init< optional <std::vector<kep_toolbox::planet::planet_ptr>, std::vector<std::vector<double> >, kep_toolbox::epoch, kep_toolbox::array3D > >())
		.def("pretty", &problem::mga_part::pretty)
		.def("get_sequence", &problem::mga_part::get_sequence)
		.add_property("vinf_in",make_function(&problem::mga_part::get_vinf_in, return_value_policy<copy_const_reference>()), &problem::mga_part::set_vinf_in,"initial incoming relative spacecraft velocity")
		.add_property("t0",make_function(&problem::mga_part::get_t0, return_value_policy<copy_const_reference>()), &problem::mga_part::set_t0, "start epoch")
		.add_property("tof",make_function(&problem::mga_part::get_tof, return_value_policy<copy_const_reference>()), &problem::mga_part::set_tof,"bounds on the times of flight for the different legs")
		.add_property("betas",&problem::mga_part::get_betas, &problem::mga_part::set_betas,"bounds on the beta angles for the different legs")
		.add_property("rps",&problem::mga_part::get_rps, &problem::mga_part::set_rps,"bounds on the periplanet heights for the different legs");

	// Travelling salesman problem, Asteroids / Debris Selection TSP (TSP-ADS)
	tsp_problem_wrapper<problem::tsp_ds>("tsp_ds","Debris Selection TSP (TSP-DS)")
		.def(init<optional<const std::vector<kep_toolbox::planet::planet_ptr>&, const std::vector<double>&, const double, const std::vector<double>&, const problem::base_tsp::encoding_type & > >())
		.def("find_subsequence", &find_subsequence_wrapper_ds)
        .add_property("planets", make_function(&problem::tsp_ds::get_planets, return_value_policy<copy_const_reference>()) )
        .add_property("values", make_function(&problem::tsp_ds::get_values, return_value_policy<copy_const_reference>()) )
        .add_property("epochs", make_function(&problem::tsp_ds::get_epochs, return_value_policy<copy_const_reference>()), "epoch schedule")
        .add_property("max_DV", &problem::tsp_ds::get_max_DV );

#ifdef PAGMO_ENABLE_GSL
	// Spheres Problems
	stochastic_problem_wrapper<problem::spheres>("mit_spheres", "Spheres problem, a neurocontroller for the MIT test-bed (absolute perception-action)")
		.def(init< optional<int,int,double,unsigned int, bool, double, std::vector<double> > >())
		.def("post_evaluate", &problem::spheres::post_evaluate)
		.def("simulate", &problem::spheres::simulate)
		.def("get_nn_weights", &problem::spheres::get_nn_weights)
		.add_property("seed",&problem::spheres::get_seed,&problem::spheres::set_seed,"Random seed used in the objective function evaluation.");

	//problem_wrapper<problem::spheres_q>("spheres_q", "Spheres problem, a neurocontroller for the MIT test-bed (body-axis perception-action)")
	//	.def(init< optional<int,int,double,unsigned int> >())
	//	.def("post_evaluate", &problem::spheres_q::post_evaluate)
	//	.def("simulate", &problem::spheres_q::simulate);
#endif
}
