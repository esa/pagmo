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

#include "../config.h"
#include "base.h"
#include "../population.h"
#include "../types.h"
#include <string>



namespace pagmo { namespace algorithm {

/// The Cross Entropy method (CE)
/**
 * The cross-entropy (CE) method attributed to Reuven Rubinstein is a general Monte Carlo approach to combinatorial and continuous multi-extremal optimization and importance sampling.
 *
 * At each call of the evolve method a number of function evaluations equal
 * to iter * pop.size() is performed.
 *
 * @see http://ie.technion.ac.il/CE/files/Misc/tutorial.pdf
 *
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 */

class __PAGMO_VISIBLE cross_entropy: public base
{
public:
	cross_entropy(int iter, double fraction_elite = 0.1);
	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
protected:
	std::string human_readable_extra() const;
private:
	const int m_iter;
	const double m_fraction_elite;
	static decision_vector	calculate_mean(std::vector<decision_vector>);
	static decision_vector calculate_std(std::vector<decision_vector>, decision_vector);
};

}} //namespaces

#endif // CROSS_ENTROY_H
