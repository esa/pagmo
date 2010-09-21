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

#ifndef PAGMO_PYTHON_BASE_ISLAND_H
#define PAGMO_PYTHON_BASE_ISLAND_H

#include <boost/python/class.hpp>
#include <stdexcept>
#include <string>

#include "../../src/algorithm/base.h"
#include "../../src/base_island.h"
#include "../../src/config.h"
#include "../../src/exceptions.h"
#include "../../src/migration/base_r_policy.h"
#include "../../src/migration/base_s_policy.h"
#include "../../src/migration/best_s_policy.h"
#include "../../src/migration/fair_r_policy.h"
#include "../../src/population.h"
#include "../../src/problem/base.h"
#include "../../src/python_locks.h"
#include "../../src/serialization.h"
#include "../utils.h"

// Forward declarations.
namespace pagmo {

class python_base_island;

}

namespace boost { namespace serialization {

template <class Archive>
void save_construct_data(Archive &, const pagmo::python_base_island *, const unsigned int);

template <class Archive>
inline void load_construct_data(Archive &, pagmo::python_base_island *, const unsigned int);

}}

namespace pagmo {

class __PAGMO_VISIBLE python_base_island:  public base_island, public boost::python::wrapper<base_island>
{
	public:
		explicit python_base_island(const problem::base &prob, const algorithm::base &algo, int n = 0,
			const double &migr_prob = 1,
			const migration::base_s_policy &s_policy = migration::best_s_policy(),
			const migration::base_r_policy &r_policy = migration::fair_r_policy()):
			base_island(prob,algo,n,migr_prob,s_policy,r_policy) {}
		explicit python_base_island(const population &pop, const algorithm::base &algo,
			const double &migr_prob = 1,
			const migration::base_s_policy &s_policy = migration::best_s_policy(),
			const migration::base_r_policy &r_policy = migration::fair_r_policy()):
			base_island(pop,algo,migr_prob,s_policy,r_policy) {}
		python_base_island(const base_island &isl):base_island(isl) {}
		python_base_island &operator=(const python_base_island &other)
		{
			base_island::operator=(other);
			return *this;
		}
		base_island_ptr clone() const
		{
			gil_state_lock lock;
			base_island_ptr retval = this->get_override("__copy__")();
			if (!retval) {
				pagmo_throw(std::runtime_error,"island's __copy__() method returns a NULL pointer, please check the implementation");
			}
			return retval;
		}
		std::string get_name() const
		{
			gil_state_lock lock;
			if (boost::python::override f = this->get_override("get_name")) {
				return f();
			}
			return base_island::get_name();
		}
		std::string default_get_name() const
		{
			return this->base_island::get_name();
		}
		population py_perform_evolution(algorithm::base_ptr a_ptr, const population &pop) const
		{
			gil_state_lock lock;
			if (boost::python::override f = this->get_override("_perform_evolution")) {
				return f(a_ptr,pop);
			}
			pagmo_throw(not_implemented_error,"island's _perform_evolution method has not been implemented");
		}
	protected:
		// An island implemented in Python is always blocking.
		bool is_blocking_impl() const
		{
			return true;
		}
		void perform_evolution(const algorithm::base &a, population &pop) const
		{
			pop = py_perform_evolution(a.clone(),pop);
		}
	private:
		template <class Archive>
		friend void boost::serialization::save_construct_data(Archive &, const python_base_island *, const unsigned int);
		template <class Archive>
		friend void boost::serialization::load_construct_data(Archive &, python_base_island *, const unsigned int);
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_island>(*this);
			ar & boost::serialization::base_object<boost::python::wrapper<base_island> >(*this);
		}
};

}

namespace boost { namespace serialization {

template <class Archive>
inline void save_construct_data(Archive &ar, const pagmo::python_base_island *isl, const unsigned int)
{
	// Save data required to construct instance.
	pagmo::problem::base_ptr prob = isl->m_pop.problem().clone();
	pagmo::algorithm::base_ptr algo = isl->m_algo->clone();
	ar << prob;
	ar << algo;
}

template <class Archive>
inline void load_construct_data(Archive &ar, pagmo::python_base_island *isl, const unsigned int)
{
	// Retrieve data from archive required to construct new instance.
	pagmo::problem::base_ptr prob;
	pagmo::algorithm::base_ptr algo;
	ar >> prob;
	ar >> algo;
	// Invoke inplace constructor to initialize instance of the island.
	::new(isl)pagmo::python_base_island(*prob,*algo);
}

}} //namespaces

BOOST_CLASS_EXPORT(pagmo::python_base_island);

#endif
