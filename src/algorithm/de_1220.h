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

#ifndef PAGMO_ALGORITHM_DE_1220_H
#define PAGMO_ALGORITHM_DE_1220_H

#include <string>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace algorithm {

/// Differential Evolution Algorithm - 1220 (our version!!)
/**
 *
 * Since Differential Evolution has always been one of PaGMO's best algorithm we
 * dared to propose our own algoritmic variant we call DE 1220. Our variant has self-adaptation of CR, F
 * and the mutation variant, so that the only parameter left to be specified is the population size.
 *
 * NOTE: when called on mixed-integer problems DE 1220 treats the integer part as fixed and optimizes
 * the continuous part.
 *
 * NOTE2: when called on stochastic optimization problems, DE 1220 changes the seed
 * at the end of each generation.
 *
 * NOTE3: the pagmo::population::individual_type::cur_v is also updated in DE 1220 as soon as a new chromosome is accepted.
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 */



class __PAGMO_VISIBLE de_1220: public base
{
public:
	de_1220(int = 100, int = 1, const std::vector<int>& = construct_default_strategies(), bool = true, double = 1e-6, double = 1e-6);
	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
private:
	static const std::vector<int> construct_default_strategies() {
			const int tmp[8] = {2,3,7,10,13,14,15,16};
			std::vector<int> retval(tmp,tmp+8);
			return retval;
	}

protected:
	std::string human_readable_extra() const;
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<int &>(m_gen);
		ar & const_cast<int &>(m_variant_adptv);
		ar & const_cast<std::vector<int> &>(m_allowed_variants);
		ar & const_cast<bool &>(m_memory);
		ar & const_cast<double &>(m_ftol);
		ar & const_cast<double &>(m_xtol);
		ar & m_f;
		ar & m_cr;
		ar & m_variants;
	}
	
	// Number of generations.
	const int m_gen;

	//Type of self-adaptive sheme for CR and F
	const int m_variant_adptv;

	//Variants allowed
	const std::vector<int> m_allowed_variants;

	//Resart option
	const bool m_memory;

	//Tolerances
	const double m_ftol;
	const double m_xtol;
	
	// Weighting factor
	mutable std::vector<double> m_f;
	
	// Crossover probability
	mutable std::vector<double> m_cr;

	// Variants of the mutation type
	mutable std::vector<int> m_variants;
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::de_1220)

#endif // DE_1220_H
