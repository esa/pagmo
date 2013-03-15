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

#ifndef PAGMO_PROBLEM_BASE_STOCHASTIC_H
#define PAGMO_PROBLEM_BASE_STOCHASTIC_H

#include "base.h"
#include "../serialization.h"
#include "../rng.h"

namespace pagmo{ namespace problem {

/// Base Stochastic Optimization Problem.
/**
 * A stochastic optimization problem is a problem that has an objective function
 * \f$ J(\mathbf x) = E_s(\mathbf x,s) \f$
 * being the average over some stochastic variable(s) \f$ s \f$. These types of problems need to
 * iherit from this base class that provides a-pseudo random number generator and a seed
 * that must be used in the objective function definition to generate pseudo random instances of
 * \f$ s\f$. These are then used to approximate the expected value with the average over n instances of s.
 * Look at pagmo::problem::inventory and pagmo::problem::spheres for typical examples.
 *
 * Optimization techniques that want to deal with these types of problems need to take care
 * of changing appropriately the seed along the optimization process as to avoid overfitting (that is
 * to avoid solving the problem only for one pseudo random sequence of s, and not for any). This must be done by
 * by a call to the change_seed() method.
 * See pagmo::algorithm::pso_stochastic
 * for a good example of such techniques.
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 */

class __PAGMO_VISIBLE base_stochastic : public base
{
	public:
		base_stochastic(int, unsigned int = 0u);
		unsigned int get_seed() const;
		void set_seed(unsigned int) const; //This is marked const as m_seed is mutable (needs to be)
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & m_drng;
			ar & m_seed;
		}

	protected:
		mutable rng_double				m_drng;
		mutable unsigned int			m_seed;
};

}} //namespaces

BOOST_SERIALIZATION_ASSUME_ABSTRACT(pagmo::problem::base_stochastic);

#endif // PAGMO_PROBLEM_BASE_STOCHASTIC_H
