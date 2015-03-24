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

#ifndef PAGMO_ALGORITHM_SA_CORANA_H
#define PAGMO_ALGORITHM_SA_CORANA_H

#include "../config.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace algorithm {

/// Simulated Annealing, Corana's version with adaptive neighbourhood.
/**
 *
 *
 * This version of A. CORANA, M. MARCHESI, C. MARTINI, and S. RIDELLA of the simulated annealing algorithm
 * is essentially an iterative random search procedure with adaptive moves along the
 * coordinate directions. It permits uphill moves under the control of metropolis criterion,
 * in the hope to avoid the first local minima encountered.
 *
 *
 * The implementation provided for PaGMO has been developed from scratch and subsequent call to the algorithm
 * are equivalent to a reannealing procedure.
 *
 * This algorithm is suitable for box-constrained single-objective continuous optimization.
 *
 * At each call of the evolve method the number of function evaluations is guaranteed to be less
 * than the total iterations as if a point is produced out of the bounds the iteration is skipped
 *
 * @see http://amcg.ese.ic.ac.uk/~jgomes/lasme/SA-corana.pdf for the original paper
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 */

class __PAGMO_VISIBLE sa_corana: public base
{
public:
	sa_corana(int niter = 1, const double &Ts = 10, const double &Tf = .1, int m_step_adj = 1, int m_bin_size = 20, const double &range = 1);
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
		ar & const_cast<int &>(m_niter);
		ar & const_cast<double &>(m_Ts);
		ar & const_cast<double &>(m_Tf);
		ar & const_cast<int &>(m_step_adj);
		ar & const_cast<int &>(m_bin_size);
		ar & const_cast<double &>(m_range);
	}  
	// Number of iterations.
	const int m_niter;
	// Starting temperature
	const double m_Ts;
	// Final temperature
	const double m_Tf;
	// Ratio of neighbourhood adjustments over temperature adjustments
	const int m_step_adj;
	// Size of the bin to evaluate the acceptance rate
	const int m_bin_size;
	// Starting neighbourhood size
	const double m_range;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::sa_corana)

#endif // PAGMO_ALGORITHM_SA_CORANA_H
