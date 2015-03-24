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

#ifndef PAGMO_PROBLEM_INVENTORY_H
#define PAGMO_PROBLEM_INVENTORY_H

#include <cstddef>
#include <string>
#include <vector>

#include "../config.h"
#include "../population.h"
#include "../rng.h"
#include "../serialization.h"
#include "../types.h"
#include "base_stochastic.h"

namespace pagmo
{
namespace problem {

/// Stochastic Programming Test Problem: Inventory Model
/**
 * This problem is a generalization of the simple inventory problem so-called of the "news-vendor",
 * widely used to introduce the main tools and techniques of stochastic programming in general
 * Assume you are a newsvendor and each week, for the next \f$ N\f$ weeks, you need to decide how many
 * journals to order (indicated with the decision variable  \f$ x_i \f$). The weekly journal demand is
 * unknown to you and is indicated with the variable \f$d_i\f$. The cost of
 * ordering journals before the week starts is \f$ c\f$, the cost of ordering journals during the week
 * (in order to meet an unforeseen demand) is \f$ b \f$ and the cost of having to hold unsold journals
 * is \f$ h \f$. The inventory level of journals will be defined by the succession:
 * \f[
 *	I_i = [I_{i-1} + x_i - d_i]_+, I_1 = 0
 * \f]
 * while the total cost of running the journal sales for \f$N\f$ weeks will be:
 * \f[
 *	J(\mathbf x, \mathbf d) = c \sum_{i=1}^N x_i+ b \sum_{i=1}^N [d_i - I_i - x_i]_+ + h \sum_{i=1}^N [I_i + x_i - d_i]_+
 * \f]
 *
 * @see www2.isye.gatech.edu/people/faculty/Alex_Shapiro/SPbook.pdf
 *
 * @author Dario Izzo (dario.izzo@esa.int)
 */
class __PAGMO_VISIBLE inventory: public base_stochastic
{
	public:
		/// Constructor from weeks, sample size and random seed
		/**
		 * Given the numer of weeks (i.e. prolem dimension), the sample size to
		 * approximate the expected value and a starting random seed, we contruct
		 * the inventory prolem
		 *
		 * @param[in] weeks integer dimension of the problem corresponding to the numer of weeks
		 * to plan the inventory for.
		 * @param[in] sample_size integer dimension of the sample used to approximate the expected value
		 * @param[in] seed unsigned integer used as starting random seed to build the pseudorandom sequences used to
		 * generate the sample
		 * @see problem::base constructors.
		 */
		inventory(int weeks = 4,int sample_size = 10, unsigned int seed = 0);
    		std::string get_name() const;
		base_ptr clone() const;
	protected:
		std::string human_readable_extra() const;
		void objfun_impl(fitness_vector &, const decision_vector &) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_stochastic>(*this);
			ar & m_weeks;
			ar & m_sample_size;
		}
	private:
		int				m_weeks;
		std::size_t			m_sample_size;

};

}
}

BOOST_CLASS_EXPORT_KEY(pagmo::problem::inventory)

#endif //PAGMO_PROBLEM_INVENTORY_H
