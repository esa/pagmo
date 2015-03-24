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

#ifndef PAGMO_UTIL_RACING_H
#define PAGMO_UTIL_RACING_H

#include <iostream>
#include <string>
#include <vector>

#include "../config.h"
#include "../serialization.h"
#include "../population.h"

///Doxygen will ignore whatever is in //! @cond
//! @cond

namespace pagmo{ namespace util {

namespace racing{
		
	class __PAGMO_VISIBLE racing_population : public population
	{
	public:
		racing_population(const population &);
		racing_population(const problem::base &);
		void set_x_noeval(const size_type, const decision_vector &);
		void set_fc(const size_type, const fitness_vector &, const constraint_vector &);
		void push_back_noeval(const decision_vector &);
		std::vector<double> get_rankings() const;
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
			stat_test_result(unsigned int N = 0): trivial(true), is_better(N, std::vector<bool>(N, false)) { }
			bool trivial;
			std::vector<std::vector<bool> > is_better;
	};

	// F-Race routines
	stat_test_result friedman_test(std::vector<racer_type> &,
	                               const std::vector<population::size_type> &,
	                               const racing_population&,
	                               double);

	stat_test_result core_friedman_test(const std::vector<std::vector<double> > &,
	                                    double delta);

	void f_race_assign_ranks(std::vector<racer_type> &,
	                         const racing_population &);

	void f_race_adjust_ranks(std::vector<racer_type> &,
	                         const std::vector<population::size_type> &);

	// Wilcoxon rank-sum routines
	stat_test_result wilcoxon_ranksum_test(std::vector<racer_type> &,
	                                       const std::vector<population::size_type> &,
	                                       const racing_population&,
	                                       double);

	stat_test_result core_wilcoxon_ranksum_test(const std::vector<std::vector<double> > &X,
	                                            double delta);
}}} //Namespaces

//! @endcond

#endif
