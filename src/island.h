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

#ifndef PAGMO_ISLAND_H
#define PAGMO_ISLAND_H

#include <string>

#include "base_island.h"
#include "config.h"
#include "algorithm/base.h"
#include "migration/base_r_policy.h"
#include "migration/base_s_policy.h"
#include "migration/best_s_policy.h"
#include "migration/fair_r_policy.h"
#include "population.h"
#include "problem/base.h"
#include "serialization.h"

// Forward declarations.
namespace pagmo {

class island;

}

namespace boost { namespace serialization {

template <class Archive>
void save_construct_data(Archive &, const pagmo::island *, const unsigned int);

template <class Archive>
inline void load_construct_data(Archive &, pagmo::island *, const unsigned int);

}}

namespace pagmo
{

/// Local island class.
/**
 * This island class will launch evolutions using local threads.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 * @author Marek Ruciński (marek.rucinski@gmail.com)
 */
class __PAGMO_VISIBLE island: public base_island
{
	public:
		island(const island &);
		explicit island(const algorithm::base &, const problem::base &, int = 0,
			const migration::base_s_policy & = migration::best_s_policy(),
			const migration::base_r_policy & = migration::fair_r_policy());
		explicit island(const algorithm::base &, const population &,
			const migration::base_s_policy & = migration::best_s_policy(),
			const migration::base_r_policy & = migration::fair_r_policy());
		island &operator=(const island &);
		base_island_ptr clone() const;
	protected:
		/** @name Evolution.
		 * Methods related to island evolution.
		 */
		//@{
		void perform_evolution(const algorithm::base &, population &) const;
		//@}
	public:
		/** @name Input/output.*/
		//@{
		std::string get_name() const;
		//@}
	private:
		template <class Archive>
		friend void boost::serialization::save_construct_data(Archive &, const island *, const unsigned int);
		template <class Archive>
		friend void boost::serialization::load_construct_data(Archive &, island *, const unsigned int);
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			// Join will be done here already.
			ar & boost::serialization::base_object<base_island>(*this);
		}
};

}

namespace boost { namespace serialization {

template <class Archive>
inline void save_construct_data(Archive &ar, const pagmo::island *isl, const unsigned int)
{
	// Save data required to construct instance.
	pagmo::algorithm::base_ptr algo = isl->m_algo->clone();
	pagmo::problem::base_ptr prob = isl->m_pop.problem().clone();
	ar << algo;
	ar << prob;
}

template <class Archive>
inline void load_construct_data(Archive &ar, pagmo::island *isl, const unsigned int)
{
	// Retrieve data from archive required to construct new instance.
	pagmo::algorithm::base_ptr algo;
	pagmo::problem::base_ptr prob;
	ar >> algo;
	ar >> prob;
	// Invoke inplace constructor to initialize instance of the algorithm.
	::new(isl)pagmo::island(*algo,*prob);
}

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::island)

#endif
