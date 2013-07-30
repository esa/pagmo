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

#ifndef PAGMO_PROBLEM_ZDT5_H
#define PAGMO_PROBLEM_ZDT5_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "base_unc_mo.h"

namespace pagmo{ namespace problem {

/// ZDT5 problem
/*
 *     This is a box-constrained integer n-dimension multi-objecive problem.
 *
 *       \f[ 
 *           F_1\left(x\right) = 1 + u \left(x_{1} \right) 
 *       \f]
 *       \f[
 *           g\left(x\right) = \sum_{i=2}^{11} v \left(u \left(x_{i} \right) \right) 
 *       \f]
 *       \f[
 *           v\left(u\left(x_{i}\right)\right) =  2 + u \left(x_{i} \right)    if u \left(x_{i} \right) < 5
 *           v\left(u\left(x_{i}\right)\right) =  1                            if u \left(x_{i} \right) = 5 
 *       \f]
 *       \f[
 *           F_2 = g \left(x \right) * 1/F_1 \left(x \right) 
 *       \f] 
 * 
 * @See = http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.30.5848&rep=rep1&type=pdf
 * @author Jainit Purohit (mjainit@gmail.com)
 */

class __PAGMO_VISIBLE zdt5 : public base_unc_mo
{
	public:
		zdt5(int = 11);
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

BOOST_CLASS_EXPORT_KEY(pagmo::problem::zdt5);

#endif // PAGMO_PROBLEM_ZDT5_H
