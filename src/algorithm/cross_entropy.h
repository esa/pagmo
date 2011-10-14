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

#ifndef PAGMO_ALGORITHM_CROSS_ENTROPY_H
#define PAGMO_ALGORITHM_CROSS_ENTROPY_H

#include <string>


#include "../config.h"
#include "base.h"
#include "../population.h"
#include "../types.h"
#include "../serialization.h"




namespace pagmo { namespace algorithm {

/// A Cross Entropy Study (CE)
/**
 * The cross-entropy (CE) method attributed to Reuven Rubinstein is a general Monte Carlo
 * approach to combinatorial and continuous multi-extremal optimization and importance sampling.
 *
 * It is similar to CMA-ES in many ways and the version here implemented is original with PaGMO developer
 * and can be considered a 'study' on sampling methods
 *
 * NOTE1: when called on mixed-integer problems CE treats the integer part as fixed and optimizes
 * the continuous part.
 *
 * NOTE2: at each call of the evolve method a number of function evaluations equal
 * to iter * pop.size() is performed.
 *
 * NOTE3: The Covariance Matrix is evaluated considering the average vector mu of the previous iteration
 * (not the current one) importing one of the key ideas of CMA-ES into this algorithm
 *
 * @see http://ie.technion.ac.il/CE/files/Misc/tutorial.pdf
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 */

class __PAGMO_VISIBLE cross_entropy: public base
{
public:
	cross_entropy(int gen = 500, double elite = 0.5, double scale = 0.05, bool screen_output = false);
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
		ar & const_cast<std::size_t &>(m_gen);
		ar & const_cast<double &>(m_elite);
		ar & const_cast<double &>(m_scale);
		ar & const_cast<bool &>(m_screen_output);
	}
	const std::size_t m_gen;
	const double m_elite;
	const double m_scale;
	const bool m_screen_output;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::cross_entropy);

#endif // CROSS_ENTROPY_H
