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

#ifndef PAGMO_PROBLEM_ZDT4_H
#define PAGMO_PROBLEM_ZDT4_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "base_unc_mo.h"

namespace pagmo{ namespace problem {

/// ZDT4 problem
/**
 *
 * This is a box-constrained continuous n-dimension multi-objecive problem.
 * \f[
 *	g\left(x\right) = 91 + \sum_{i=2}^{n} \left[x_i^2 - 10 \cos \left(4\pi x_i \right) \right]
 * \f]
 * \f[
 * 	F_1 \left(x\right) = x_1
 * \f]
 * \f[
 *      F_2 \left(x\right) = g(x) \left[ 1 - \sqrt{x_1 / g(x)} \right]  x_1 \in [0,1], x_i \in \left[ -5,5 \right] i=2, \cdots, 10.
 *
 * \f]
 *
 * @see http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.18.4257&rep=rep1&type=pdf
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 */

class __PAGMO_VISIBLE zdt4 : public base_unc_mo
{
	public:
		zdt4(size_type = 10);
		base_ptr clone() const;
		std::string get_name() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		double convergence_metric(const decision_vector &) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_unc_mo>(*this);
		}
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::zdt4);

#endif // PAGMO_PROBLEM_ZDT4_H
