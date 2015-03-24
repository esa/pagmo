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

#ifndef PAGMO_PROBLEM_LAPLACE_H
#define PAGMO_PROBLEM_LAPLACE_H

#include <string>
#include <vector>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "../AstroToolbox/mga_dsm.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// Laplace problem
/**
 * \image html laplace.jpg "Laplace: a jupiter mission"
 * \image latex laplace.jpg "Laplace: a jupiter mission" width=5cm
 *
 * Laplace is a Jupiter mission proposed as a candidate mission for the Cosmic Vision 2015-2025
 * and has been one of the mission concepts selected by the SSAC in October 2007 for an Assessment
 * study. The SSAC had decided that an Outer Planet mission should be one of the L-class
 * mission candidates for the first slice of the Cosmic Vision plan, and recommended that
 * both Laplace and Tandem concepts should be studied initially, with a decision on which
 * one to pursue to take place following a better definition of the mission's characteristics.
 *
 * In early February 2009 ESA and NASA jointly announced that the Europa-Jupiter
 * System Mission, or Laplace, would be the candidate for the first L mission,
 * with a possible launch date in 2020.
 *
 * We transcribe here Laplace interplanetary trajectory transfer problem as an MGA-DSM
 * problem leaving to the user the possibility to specify the fly-by sequence.
 * This creates a problem of dimension \f$6 + 4(n-1)\f$, where \f$n\f$ is the number of
 * trajectory legs. Objective function is the total DV with a 200m/s penalty per month past the 8yr
 * of flight time.
 *
 *
 * @see http://sci.esa.int/science-e/www/area/index.cfm?fareaid=107
 * @author Dario Izzo (dario.izzo@esa.int)
 */
class __PAGMO_VISIBLE laplace: public base
{
		static const int default_sequence[5];
		static const int* get_default_sequence();
	public:
		laplace(const std::vector<int> & = std::vector<int>(get_default_sequence(),get_default_sequence() + 5));
		base_ptr clone() const;
		std::string get_name() const;
		std::string pretty(const std::vector<double> &x) const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		//void set_sparsity(int &, std::vector<int> &, std::vector<int> &) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & problem;
		}
		mgadsmproblem problem;
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::problem::laplace)

#endif // PAGMO_PROBLEM_LAPLACE_H
