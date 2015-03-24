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

#ifndef PAGMO_ALGORITHM_JDE_H
#define PAGMO_ALGORITHM_JDE_H

#include <string>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace algorithm {

/// jDE - Differential Evolution Algorithm - Self-Adaptive C and R (2011)
/**
 *
 * \image html de.jpg "Differential Evolution block diagram."
 * \image latex de.jpg "Differential Evolution block diagram." width=5cm
 *
 * Since its creation, the original Differential Evolution (pagmo::algorithm::de) algorithm
 * has been modified several times and many improvements have been suggested. We thus provide in PaGMO, together with the
 * original version of the algorithm, a modern version of the algorithm, with self-adaptation of its 
 * parameters pagmo::algorithm::de::m_cr and pagmo::algorithm::de::m_f and some of more recombination variants.
 *
 * NOTE: when called on mixed-integer problems DE treats the integer part as fixed and optimizes
 * the continuous part.
 *
 * NOTE2: when called on stochastic optimization problems, DE changes the seed
 * at the end of each generation.
 *
 * NOTE3: the pagmo::population::individual_type::cur_v is also updated along DE as soon as a new chromosome is accepted.
 *
 *
 * @see http://labraj.uni-mb.si/images/0/05/CEC09_slides_Brest.pdf  where m_variant_adptv = 1 is studied.
 * @see http://sci2s.ugr.es/EAMHCO/pdfs/contributionsCEC11/05949732.pdf for a paper where a similar apporach to m_variant_adptv=2 is described
 * 'modern' de version are used.
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 */

class __PAGMO_VISIBLE jde : public base
{
public:
	jde(int = 100, int = 2, int = 1, double = 1e-6, double = 1e-6, bool = false);
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
		ar & const_cast<int &>(m_gen);
		ar & const_cast<int &>(m_variant);
		ar & const_cast<int &>(m_variant_adptv);
		ar & const_cast<double &>(m_ftol);
		ar & const_cast<double &>(m_xtol);
		ar & m_f;
		ar & m_cr;
		ar & const_cast<bool &>(m_memory);
	}
	
	// Number of generations.
	const int m_gen;
	
	// Weighting factor
	mutable std::vector<double> m_f;
	
	// Crossover probability
	mutable std::vector<double> m_cr;
	
	// Algoritmic variant
	const int m_variant;
	// Self-adaptation strategy
	const int m_variant_adptv;
	const double m_ftol;
	const double m_xtol;

	// Memory option
	const bool m_memory;
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::jde)

#endif // JDE_H
