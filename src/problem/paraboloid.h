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

// 01/02/10 Created by Francesco Biscani.

#ifndef PAGMO_PROBLEM_PARABOLOID_H
#define PAGMO_PROBLEM_PARABOLOID_H

#include <cstddef>
#include <string>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"

namespace pagmo { namespace problem {

/// N-dimensional paraboloid problem.
/**
 * \image html paraboloid.png "Two-dimensional paraboloid."
 * \image latex paraboloid.png "Two-dimensional paraboloid." width=5cm
 *
 * This is a box-constrained continuous single-objecive problem, also known as "sphere function", "sphere model" and "De Jong's function 1". 
 * The objective function for an N-dimensional instance of this problem is:
 * \f[
 * 	f\left(x_1,\ldots,x_N \right) = \sum_{i=1}^N x_i^2.
 * \f]
 * The only and global minimum is in the origin. Useful for testing/benchmarking purposes.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE paraboloid: public base
{
	public:
		paraboloid();
		/// Constructor from raw arrays.
		/**
		 * @see problem::base constructors.
		 */
		template <std::size_t N>
		paraboloid(const double (&v1)[N], const double (&v2)[N]):base(v1,v2) {}
		/// Constructor from iterators.
		/**
		 * @see problem::base constructors.
		 */
		template <class Iterator1, class Iterator2>
		paraboloid(Iterator1 start1, Iterator1 end1, Iterator2 start2, Iterator2 end2):base(start1,end1,start2,end2) {}
		paraboloid(const decision_vector &, const decision_vector &);
		base_ptr clone() const;
		std::string get_name() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
		}
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::problem::paraboloid);

#endif
