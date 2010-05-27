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
#include <boost/python/pure_virtual.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <string>

#include "../../src/algorithms.h"
#include "../../src/population.h"
#include "../exceptions.h"

using namespace boost::python;
using namespace pagmo;

// Wrapper method for algorithm evolve that uses copy instead of pass-by-non-const-reference.
static inline population evolve_copy(const algorithm::base &a, const population &pop)
{
	population pop_copy(pop);
	a.evolve(pop_copy);
	return pop_copy;
}

template <class Algorithm>
static inline class_<Algorithm,bases<algorithm::base> > algorithm_wrapper(const char *name, const char *descr)
{
	class_<Algorithm,bases<algorithm::base> > retval(name,descr,init<const Algorithm &>());
	retval.def("__copy__", &Algorithm::clone);
	retval.def("evolve", &evolve_copy);
	return retval;
}

struct python_algorithm: algorithm::base, wrapper<algorithm::base>
{
	python_algorithm():algorithm::base() {}
	python_algorithm(const algorithm::base &p):algorithm::base(p) {}
	algorithm::base_ptr clone() const
	{
		return this->get_override("__copy__")();
	}
	void evolve(population &p) const
	{
		p = py_evolve(p);
	}
	population py_evolve(const population &p) const
	{
		return this->get_override("evolve")(p);
	}
	bool is_blocking() const
	{
		return true;
	}
	std::string human_readable_extra() const
	{
		return py_human_readable_extra();
	}
	std::string py_human_readable_extra() const
	{
		if (override f = this->get_override("_human_readable_extra")) {
			return f();
		}
		return algorithm::base::human_readable_extra();
	}
	std::string get_name() const
	{
		if (override f = this->get_override("get_name")) {
			return f();
		}
		return algorithm::base::get_name();
	}
	std::string default_get_name() const
	{
		return this->algorithm::base::get_name();
	}
};

BOOST_PYTHON_MODULE(_algorithm) {
	// Translate exceptions for this module.
	translate_exceptions();

	// Expose base algorithm class, including the virtual methods.
	class_<python_algorithm>("_base",init<>())
		.def(init<const algorithm::base &>())
		.def("__repr__", &algorithm::base::human_readable)
		.def("is_blocking",&algorithm::base::is_blocking)
		// Virtual methods that can be (re)implemented.
		.def("__copy__",pure_virtual(&algorithm::base::clone))
		.def("get_name",&algorithm::base::get_name,&python_algorithm::default_get_name)
		.def("evolve",&python_algorithm::py_evolve)
		.def("_human_readable_extra",&python_algorithm::py_human_readable_extra);

	// Expose algorithms.

	// Null.
	algorithm_wrapper<algorithm::null>("null","Null algorithm.")
		.def(init<>());

	// IHS.
	algorithm_wrapper<algorithm::ihs>("ihs","Improved harmony search.")
		.def(init<int, optional<const double &, const double &, const double &, const double &, const double &> >());
	
	// CS.
	algorithm_wrapper<algorithm::cs>("cs","Compass search solver.")
		.def(init<const int &, const double &, optional<const double &, const double &> >());

	// Monte-carlo.
	algorithm_wrapper<algorithm::monte_carlo>("monte_carlo","Monte-Carlo search.")
		.def(init<int>());

	// Artificial Bee Colony Optimization (ABC).
	algorithm_wrapper<algorithm::bee_colony>("bee_colony","Artificial Bee Colony optimization (ABC) algorithm.")
		.def(init<int,optional<int> >());
	
	// Monotonic Basin Hopping.
	algorithm_wrapper<algorithm::mbh>("mbh","Monotonic Basin Hopping.")
		.def(init<const algorithm::base &,optional<int, double> >());
	
	// Multistart.
	algorithm_wrapper<algorithm::ms>("ms","Multistart.")
		.def(init<const algorithm::base &, int >());
	
	// Particle Swarm Optimization.
	algorithm_wrapper<algorithm::pso>("pso","Particle swarm optimization.")
		.def(init<int,optional<double, double, double, double, int> >());
	
	// Simple Genetic Algorithm.
	algorithm_wrapper<algorithm::sga>("sga","Simple Genetic Algorithm.")
		.def(init<int, const double &, const double &, optional<int, algorithm::sga::mutation::type, double, algorithm::sga::selection::type, algorithm::sga::crossover::type> >());
	
	// Differential evolution.
	algorithm_wrapper<algorithm::de>("de","Differential evolution algorithm.")
		.def(init<int,optional<const double &, const double &, int> >());
	// ASA.
	algorithm_wrapper<algorithm::sa_corana>("sa_corana","Simulated annealing, Corana's version with adaptive neighbourhood.")
		.def(init<int, const double &, const double &, optional<int,int,const double &> >());

	// GSL algorithms.
	#ifdef PAGMO_ENABLE_GSL

	// GSL's BFGS.
	algorithm_wrapper<algorithm::gsl_bfgs>("gsl_bfgs","GSL BFGS algorithm.")
		.def(init<optional<int, const double &, const double &, const double &, const double &> >());

	algorithm_wrapper<algorithm::gsl_bfgs2>("gsl_bfgs2","GSL BFGS2 algorithm.")
		.def(init<optional<int, const double &, const double &, const double &, const double &> >());
	
	// GSL's Fletcher-Reeves.
	algorithm_wrapper<algorithm::gsl_fr>("gsl_fr","GSL Fletcher-Reeves algorithm.")
		.def(init<optional<int, const double &, const double &, const double &, const double &> >());

	// GSL's Nelder-Mead.
	algorithm_wrapper<algorithm::gsl_nm>("gsl_nm","GSL Nelder-Mead simplex method.")
		.def(init<optional<int, const double &, const double &> >());

	// GSL's Nelder-Mead, version 2.
	algorithm_wrapper<algorithm::gsl_nm2>("gsl_nm2","GSL Nelder-Mead simplex method, version 2.")
		.def(init<optional<int, const double &, const double &> >());

	// GSL's Nelder-Mead, version 2 + random initial simplex.
	algorithm_wrapper<algorithm::gsl_nm2rand>("gsl_nm2rand","GSL Nelder-Mead simplex method, version 2 with randomised initial simplex.")
		.def(init<optional<int, const double &, const double &> >());

	// GSL's Polak-Ribiere.
	algorithm_wrapper<algorithm::gsl_pr>("gsl_pr","GSL Polak-Ribiere algorithm.")
		.def(init<optional<int, const double &, const double &, const double &, const double &> >());

	#endif

	// NLopt algorithms.
	#ifdef PAGMO_ENABLE_NLOPT

	// NLopt's COBYLA.
	algorithm_wrapper<algorithm::nlopt_cobyla>("nlopt_cobyla","NLopt's COBYLA algorithm.")
		.def(init<optional<int, const double &> >());

	// NLopt's BOBYQA.
	algorithm_wrapper<algorithm::nlopt_bobyqa>("nlopt_bobyqa","NLopt's BOBYQA algorithm.")
		.def(init<optional<int, const double &> >());

	// NLopt's Sbplx.
	algorithm_wrapper<algorithm::nlopt_sbplx>("nlopt_sbplx","NLopt's Sbplx algorithm.")
		.def(init<optional<int, const double &> >());

	#endif

	// Snopt solver.
	#ifdef PAGMO_ENABLE_SNOPT
	algorithm_wrapper<algorithm::snopt>("snopt","Snopt solver.")
		.def(init<int, optional<double, double> >());
	
	#endif

	// Register to_python conversion from smart pointer.
	register_ptr_to_python<algorithm::base_ptr>();
}
