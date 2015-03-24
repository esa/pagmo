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

#ifndef PAGMO_ALGORITHM_CS_H
#define PAGMO_ALGORITHM_CS_H

#include "../config.h"
#include "../problem/base.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace algorithm {

/// The Compass Search Solver (CS)
/**
 * In the review paper by Kolda, Lewis, Torczon: 'Optimization by Direct Search: New Perspectives on Some Classical and Modern Methods'
 * published in the SIAM Journal Vol. 45, No. 3, pp. 385-482 (2003) we read the following description of the compass search
 * algorithm:
 * 
 * 'Davidon describes what is one of the earliest examples of a direct
 * search method used on a digital computer to solve an optimization problem:
 * Enrico Fermi and Nicholas Metropolis used one of the first digital computers,
 * the Los Alamos Maniac, to determine which values of certain theoretical
 * parameters (phase shifts) best fit experimental data (scattering cross
 * sections). They varied one theoretical parameter at a time by steps
 * of the same magnitude, and when no such increase or decrease in any one
 * parameter further improved the fit to the experimental data, they halved
 * the step size and repeated the process until the steps were deemed sufficiently
 * small. Their simple procedure was slow but sure, and several of us
 * used it on the Avidac computer at the Argonne National Laboratory for
 * adjusting six theoretical parameters to fit the pion-proton scattering data
 * we had gathered using the University of Chicago synchrocyclotron.
 * While this basic algorithm undoubtedly predates Fermi and Metropolis, it has remained
 * a standard in the scientific computing community for exactly the reason observed
 * by Davidon: it is slow but sure'.
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 *
 * @see http://www.cs.wm.edu/~va/research/sirev.pdf
 */

class __PAGMO_VISIBLE cs: public base
{
public:
	cs(const int& max_eval = 1, const double &stop_range = 0.01,const double &start_range = 0.1, const double &reduction_coeff = 0.5 );
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
		ar & const_cast<double &>(m_stop_range);
		ar & const_cast<double &>(m_start_range);
		ar & const_cast<double &>(m_reduction_coeff);
		ar & const_cast<double &>(m_max_eval);
	}  
	// Stopping search length
	const double m_stop_range;
	// Starting search length
	const double m_start_range;
	// Search length reduction coefficient
	const double m_reduction_coeff;
	// Maximum function evaluations
	const double m_max_eval;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::cs)

#endif // PAGMO_ALGORITHM_CS_H
