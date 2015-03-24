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

#ifndef PAGMO_PROBLEM_BUKINF6_H
#define PAGMO_PROBLEM_BUKINF6_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"

namespace pagmo { namespace problem {

/// Bukin f6 function
/**
 * \image html bukin.jpg "Plot of Bukin's f6 function."
 * \image latex bukin.jpg "Plot of Bukin's f6 function." width=3cm
 *
 * Bukin’s functions are almost fractal (with fine seesaw edges) in 
 * the surroundings of their minimal points. Due to this property, they are extremely 
 * difficult to optimize by any method of global (or local) optimization and find correct 
 * values of decision variables (i.e. xi
 * for i=1,2).
 *
 * The function is defined by:
 * \f[
 * 	f(x, y) = 100\sqrt{\left|y - 0.01 x^2\right|} + 0.01\left|x+10\right|
 * \f]
 * The function has one local maximum at \f$ x = -10 \f$ and \f$ y = 1 \f$, where \f$ f(x,y) = 0 \f$
 *
 * The bounds of the search space are \f$ lb = [-15,-3], ub = [-5,3] \f$ range.
 *
 * @see http://mpra.ub.uni-muenchen.de/1743/1/MPRA_paper_1743.pdf
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 */
class __PAGMO_VISIBLE bukin: public base
{
	public:
		bukin();
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

BOOST_CLASS_EXPORT_KEY(pagmo::problem::bukin)

#endif
