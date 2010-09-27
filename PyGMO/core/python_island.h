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

#ifndef PAGMO_PYTHON_ISLAND_H
#define PAGMO_PYTHON_ISLAND_H

#include "../../src/config.h"
#include "../../src/base_island.h"
#include "../../src/island.h"
#include "../../src/algorithm/base.h"
#include "../../src/migration/base_r_policy.h"
#include "../../src/migration/base_s_policy.h"
#include "../../src/migration/best_s_policy.h"
#include "../../src/migration/fair_r_policy.h"
#include "../../src/population.h"
#include "../../src/problem/base.h"
#include "../../src/serialization.h"

// Forward declarations.
namespace pagmo {

class python_island;

}

namespace boost { namespace serialization {

template <class Archive>
void save_construct_data(Archive &, const pagmo::python_island *, const unsigned int);

template <class Archive>
inline void load_construct_data(Archive &, pagmo::python_island *, const unsigned int);

}}

namespace pagmo {

// We need to reimplement the local island class to make it conditionally blocking in Python.
class __PAGMO_VISIBLE python_island: public island
{
	public:
		explicit python_island(const problem::base &prob, const algorithm::base &algo, int n = 0,
			const double &migr_prob = 1,
			const migration::base_s_policy &s_policy = migration::best_s_policy(),
			const migration::base_r_policy &r_policy = migration::fair_r_policy()):
			island(prob,algo,n,migr_prob,s_policy,r_policy) {}
		explicit python_island(const population &pop, const algorithm::base &algo,
			const double &migr_prob = 1,
			const migration::base_s_policy &s_policy = migration::best_s_policy(),
			const migration::base_r_policy &r_policy = migration::fair_r_policy()):
			island(pop,algo,migr_prob,s_policy,r_policy) {}
		python_island &operator=(const python_island &other)
		{
			island::operator=(other);
			return *this;
		}
		base_island_ptr clone() const
		{
			return base_island_ptr(new python_island(*this));
		}
	protected:
		bool is_blocking_impl() const
		{
			return false;
		}
	private:
		template <class Archive>
		friend void boost::serialization::save_construct_data(Archive &, const python_island *, const unsigned int);
		template <class Archive>
		friend void boost::serialization::load_construct_data(Archive &, python_island *, const unsigned int);
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			// Join will be done here already.
			ar & boost::serialization::base_object<island>(*this);
		}
};

}

namespace boost { namespace serialization {

template <class Archive>
inline void save_construct_data(Archive &ar, const pagmo::python_island *isl, const unsigned int)
{
	// Save data required to construct instance.
	pagmo::problem::base_ptr prob = isl->m_pop.problem().clone();
	pagmo::algorithm::base_ptr algo = isl->m_algo->clone();
	ar << prob;
	ar << algo;
}

template <class Archive>
inline void load_construct_data(Archive &ar, pagmo::python_island *isl, const unsigned int)
{
	// Retrieve data from archive required to construct new instance.
	pagmo::problem::base_ptr prob;
	pagmo::algorithm::base_ptr algo;
	ar >> prob;
	ar >> algo;
	// Invoke inplace constructor to initialize instance of the island.
	::new(isl)pagmo::python_island(*prob,*algo);
}

}} //namespaces

BOOST_CLASS_EXPORT(pagmo::python_island);

#endif
