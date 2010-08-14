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

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/vector.hpp>
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
	public:
		mpi_island(const mpi_island &);
		explicit mpi_island(const problem::base &, const algorithm::base &, int = 0, int = 0,
			const double & = 1,
			const migration::base_s_policy & = migration::best_s_policy(),
			const migration::base_r_policy & = migration::fair_r_policy());
		explicit mpi_island(const population &, const algorithm::base &, int = 0,
			const double & = 1,
			const migration::base_s_policy & = migration::best_s_policy(),
			const migration::base_r_policy & = migration::fair_r_policy());
		mpi_island &operator=(const mpi_island &);
		base_island_ptr clone() const;
		/** @name Evolution.
		 * Methods related to island evolution.
		 */
		//@{
		bool is_thread_blocking() const;
	protected:
		void perform_evolution(const algorithm::base &, population &) const;
		//@}
	public:
		// the id of the processor that will perform the evolution of the island object
		int processor_id;
		/** @name Input/output.*/
		//@{
		std::string get_name() const;
		//@}
	private:
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive &ar, const unsigned int version){
		    std::cout << "de-/serializing mpi_island " << version << std::endl;
			ar & boost::serialization::base_object<base_island>(*this);			
		}
};

}

#endif
