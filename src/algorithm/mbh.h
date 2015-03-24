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

#ifndef PAGMO_ALGORITHM_MBH_H
#define PAGMO_ALGORITHM_MBH_H

#include <string>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"
#include "cs.h"



namespace pagmo { namespace algorithm {

/// Monotonic Basin Hopping (generalized)
/**
 *
 * \image html mbh.png "A schematic diagram illustrating the effects of the basin hopping transformation on a one dimensional landscape."
 * \image latex mbh.png "A schematic diagram illustrating the effects of the basin hopping transformation on a one dimensional landscape." width=5cm
 *
 * Monotonic basin hopping, or simply, basin hopping, is an algorithm rooted in the idea of transforming each
 * the objective function at \f$\mathbf x\f$ into the local minima found starting from \f$\mathbf x\f$
 * This simple idea allowed a substantial increase of efficiency in solving problems, such as the Lennard-Jones
 * cluster or the MGA-DSM interplanetary trajectory problem that are conjectured to have a so-called funnel structure.
 *
 * In PaGMO we provide an original generalization of the method that operates on any pagmo::population using any pagmo::algorithm.
 * When a population containing a single individual is used and coupled with a local optimizer, the original method is recovered.
 * The pseudo code of the generalized version is:
@verbatim
> Select a pagmo::population
> Select a pagmo::algorithm
> Store best individual
> while i < stop_criteria
> > Perturb the population in a selected neighbourhood
> > Evolve the population using the algorithm
> > if the best individual is improved (according to the problem::compare_fc criteria)
> > > increment i
> > > update best individual
> > else
> > > i = 0

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
	mbh(const base & = cs(), int stop = 5, double perturb = 5e-2);
	mbh(const base &, int stop, const std::vector<double> &perturb);
	mbh(const mbh &);
	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
	base_ptr get_algorithm() const;
	void set_algorithm(const base &);
protected:
	std::string human_readable_extra() const;
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & m_local;
		ar & const_cast<int &>(m_stop);
		ar & m_perturb;
	}
	base_ptr m_local;
	// Consecutive non improving iterations
	const int m_stop;
	// Perturbation of the population
	mutable std::vector<double> m_perturb;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::mbh)

#endif // PAGMO_ALGORITHM_MBH_H
