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

#ifndef PAGMO_MPI_ISLAND_H
#define PAGMO_MPI_ISLAND_H

#include <boost/scoped_ptr.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <set>
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

class mpi_island;

}

namespace boost { namespace serialization {

template <class Archive>
void save_construct_data(Archive &, const pagmo::mpi_island *, const unsigned int);

template <class Archive>
inline void load_construct_data(Archive &, pagmo::mpi_island *, const unsigned int);

}}

namespace pagmo
{

/// mpi island class.
/**
 * This island class will launch evolutions using local threads.
 *
 * @author Dante Stroe (dante.stroe@gmail.com)
 */
class __PAGMO_VISIBLE mpi_island: public base_island
{
		typedef boost::lock_guard<boost::mutex> lock_type;
		template <class Archive>
		friend void boost::serialization::save_construct_data(Archive &, const pagmo::mpi_island *, const unsigned int);
		template <class Archive>
		friend void boost::serialization::load_construct_data(Archive &, pagmo::mpi_island *, const unsigned int);
	public:
		mpi_island(const mpi_island &);
		explicit mpi_island(const problem::base &, const algorithm::base &, int = 0,
			const double & = 1,
			const migration::base_s_policy & = migration::best_s_policy(),
			const migration::base_r_policy & = migration::fair_r_policy());
		explicit mpi_island(const population &, const algorithm::base &,
			const double & = 1,
			const migration::base_s_policy & = migration::best_s_policy(),
			const migration::base_r_policy & = migration::fair_r_policy());
		mpi_island &operator=(const mpi_island &);
		base_island_ptr clone() const;
	protected:
		/** @name Evolution.
		 * Methods related to island evolution.
		 */
		//@{
		bool is_blocking_impl() const;
		void perform_evolution(const algorithm::base &, population &) const;
		//@}
	public:
		/** @name Input/output.*/
		//@{
		std::string get_name() const;
		//@}
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			// Join is already done in base_island.
			ar & boost::serialization::base_object<base_island>(*this);
		}
		static void init_processors();
		static int acquire_processor();
		static void release_processor(int);
	private:
		static boost::mutex				m_mutex;
		static boost::scoped_ptr<std::set<int> >	m_available_processors;
};

}

namespace boost { namespace serialization {

template <class Archive>
inline void save_construct_data(Archive &ar, const pagmo::mpi_island *isl, const unsigned int)
{
	// Save data required to construct instance.
	pagmo::problem::base_ptr prob = isl->m_pop.problem().clone();
	pagmo::algorithm::base_ptr algo = isl->m_algo->clone();
	ar << prob;
	ar << algo;
}

template <class Archive>
inline void load_construct_data(Archive &ar, pagmo::mpi_island *isl, const unsigned int)
{
	// Retrieve data from archive required to construct new instance.
	pagmo::problem::base_ptr prob;
	pagmo::algorithm::base_ptr algo;
	ar >> prob;
	ar >> algo;
	// Invoke inplace constructor to initialize instance of the algorithm.
	::new(isl)pagmo::mpi_island(*prob,*algo);
}

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::mpi_island);

#endif
