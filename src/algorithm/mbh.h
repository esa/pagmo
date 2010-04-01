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

#ifndef PAGMO_ALGORITHM_MBH_H
#define PAGMO_ALGORITHM_MBH_H

#include <string>

#include "../config.h"
#include "../population.h"
#include "base.h"


namespace pagmo { namespace algorithm {

/// Monotonic Basin Hopping (generalized)
/**
 *
 * \image html mbh.png "A schematic diagram illustrating the eﬀects of the basin hopping transformation on a one dimensional landscape".
 * \image latex mbh.png "A schematic diagram illustrating the eﬀects of the basin hopping transformation on a one dimensional landscape". width=5cm
 *
 * Monotonic basin hopping, or simply, basin hopping, is an algorithm rooted in the idea of transforming each
 * the objective function at \f$\mathbf x\f$ into the local minima found starting from \f$\mathbf x\f$
 * This simple idea allowed a substantial increase of efficiency in solving problems, such as the Lennard-Jones
 * cluster or the MGA-DSM interplanetary trajectory problem that are conjectured to have a so-called funnel structure.
 *
 * In PaGMO we provide an original generalization of the method to operate on pagmo::population and using any pagmo::algorithm.
 * The pseudo code of the generalized version is:
@verbatim
> Select a pagmo::population
> Select a pagmo::algorithm
> while i < stop_criteria
> > Evolve the population using the algorithm
> > Increment i if the champion is not improved
> > i = 0 if a better champion is found
> > Perturb the population in a selected neighbourhood
@endverbatim
 *
 *
 * @see http://arxiv.org/pdf/cond-mat/9803344 for the paper inroducing the basin hopping idea for a Lennard-Jones cluster optimization
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 */
		
class __PAGMO_VISIBLE mbh: public base
{
public:
	mbh(base local, int stop = 50, double perturb = 5e-2);
	base_ptr clone() const;
	void evolve(population &) const;
protected:
	std::string human_readable_extra() const;
private:
	base m_local;
	// Consecutive non improving iterations
	const int m_stop;
	// Perturbation of the population
	const double m_perturb;

};

}} //namespaces

#endif // PAGMO_ALGORITHM_MBH_H
