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
#include "base.h"

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
class __PAGMO_VISIBLE inventory: public base
{
	public:
		inventory(int = 4,int = 10);
		base_ptr clone() const;
	protected:
		bool equality_operator_extra(const base &) const;
		void pre_evolution(population &) const;
		void objfun_impl(fitness_vector &, const decision_vector &) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & m_seed;
			ar & m_weeks;
			ar & m_sample_size;
			ar & m_drng;
		}
	private: // save_construct_data needs to be able to access these attributes
		mutable int			m_seed;
	public:
		int				m_weeks;
		std::size_t			m_sample_size;
	private:
		mutable rng_double		m_drng;
};

}
}

BOOST_CLASS_EXPORT_KEY(pagmo::problem::inventory);

#endif //PAGMO_PROBLEM_INVENTORY_H
