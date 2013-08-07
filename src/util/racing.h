/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
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

#ifndef PAGMO_UTIL_RACING_H
#define PAGMO_UTIL_RACING_H

#include <iostream>
#include <string>
#include <vector>

#include "../config.h"
#include "../serialization.h"
#include "../population.h"

namespace pagmo{ namespace util {

namespace racing{
	
///Doxygen will ignore whatever is in //! @cond
//! @cond
	
	class __PAGMO_VISIBLE racing_population : public population
	{
	public:
		racing_population(const population &);
		void set_x_noeval(const size_type, const decision_vector &);
		void set_fc(const size_type, const fitness_vector &, const constraint_vector &);
	};

	struct racer_type
	{
		public:
			racer_type(): m_mean(0), active(false) { }

			// Using double type to cater for tied ranks
			std::vector<double> m_hist;
			double m_mean;
			bool active;

			unsigned int length()
			{
				return m_hist.size();
			}

			void reset()
			{
				m_hist.clear();
				m_mean = 0;
				active = false;
			}


		private:
			friend class boost::serialization::access;
			template <class Archive>
			void serialize(Archive &ar, const unsigned int)
			{
				ar & m_hist;
				ar & m_mean;
				ar & active;
			}
	};

	struct stat_test_result{
		public:
			stat_test_result(): trivial(true), is_better(0, std::vector<bool>(0, false)) { }
			bool trivial;
			std::vector<std::vector<bool> > is_better;
	};

	// TODO: May create a base class statistical_test_base and then
	// let specific stat tests derive from there.
	
	// Used in F-Race
	stat_test_result friedman_test(const std::vector<std::vector<double> > &, double delta);

	// Used in Hoeffding / Bernstein race
	enum bound_type {HOEFFDING = 0, BERNSTEIN = 1};
	stat_test_result bound_based_test(const std::vector<std::vector<double> > &, double delta, bound_type);

	void f_race_assign_ranks(std::vector<racer_type> &, const racing_population &);
	void f_race_adjust_ranks(std::vector<racer_type> &, const std::vector<population::size_type> &);

//! @endcond
}

}}

#endif
