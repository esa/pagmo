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

#ifndef PAGMO_ALGORITHM_FIREFLY_H
#define PAGMO_ALGORITHM_FIREFLY_H

#include <string>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace algorithm {

/// The Firefly algorithm
/**
 * The firefly algorithm (FA) is a metaheuristic algorithm, inspired by the flashing behaviour of fireflies.
 *
 * At each call of the evolve method a number of function evaluations equal
 * to gen * pop.size() * pop.size() is performed.
 *
 * NOTE: when called on mixed-integer problems Firefly treats the integer part as fixed and optimizes
 * the continuous part.
 *
 * The algorithm used is the one provided in the first paper attached with 2 differences:
 *
 * \f$ \gamma \f$ is calculated with the adaptive method explained in the second paper attached (Equation 3)
 * The attactiveness is calculated as \f$ \beta = \beta_0 \exp(-\gamma r) \f$
 *
 * These modifications make the algorithm perform better on all the test problems we experimented.
 *
 * @see http://arxiv.org/abs/1003.1466
 * @see http://www.springerlink.com/content/au3w21311g465007/
 *
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 */

class __PAGMO_VISIBLE firefly: public base
{
public:
	firefly(int gen = 1, double alpha = 0.01, double beta = 1.0, double gamma = 0.8);
	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
protected:
	std::string human_readable_extra() const;
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<int &>(m_iter);
		ar & const_cast<double &>(m_alpha);
		ar & const_cast<double &>(m_beta);
		ar & const_cast<double &>(m_gamma);
	}
	const int m_iter;
	const double m_alpha;
	const double m_beta;
	const double m_gamma;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::firefly)

#endif // PAGMO_ALGORITHM_FIREFLY_H
