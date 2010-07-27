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

// 04/01/2009: Initial version by Francesco Biscani.

#ifndef PAGMO_ISLAND_H
#define PAGMO_ISLAND_H

#include <boost/scoped_ptr.hpp>
#include <boost/thread/thread.hpp>
#include <iostream>
#include <string>

#include "base_island.h"
#include "config.h"
#include "algorithm/base.h"
#include "migration/base_r_policy.h"
#include "migration/base_s_policy.h"
#include "migration/best_s_policy.h"
#include "migration/fair_r_policy.h"
#include "problem/base.h"

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
		explicit island(const problem::base &, const algorithm::base &, int = 0,
			const double & = 1,
			const migration::base_s_policy & = migration::best_s_policy(),
			const migration::base_r_policy & = migration::fair_r_policy());
		explicit island(const population &, const algorithm::base &,
			const double & = 1,
			const migration::base_s_policy & = migration::best_s_policy(),
			const migration::base_r_policy & = migration::fair_r_policy());
		island &operator=(const island &);
		~island();
		/** @name Evolution.
		 * Methods related to island evolution.
		 */
		//@{
		bool is_blocking() const;
		//@}
		/** @name Input/output.*/
		//@{
		std::string get_name() const;
		//@}
};

}

#endif
