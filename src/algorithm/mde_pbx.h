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

#ifndef PAGMO_ALGORITHM_MDE_PBX_H
#define PAGMO_ALGORITHM_MDE_PBX_H

#include <string>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace algorithm {

/// MDE_pBX - Differential Evolution variant
/**
 *
 * MDE_pBX is a member of the Differential Evolution (pagmo::algorithm::de) algorithms.
 * It uses the mutation scheme DE/current-to-gr_best/1 which takes the best individual on a q%
 * sample of the complete population for the generation of the donor vector (by default: q% = 15%). 
 * The crossover scheme used is p-best crossover which uses a random vector choosen uniformly 
 * out of the p best individual from the population. The fitness scale Factor F and 
 * Crossover Probability Cr are adapted each generation and thus do not need to be provided as input.
 * For each individual, Cr and F are drawn from a random distributions, defined by Cr_m and F_m.
 * Each time, an individual with better fitness is created, the corresponding values for F and Cr
 * are memorized. After one evolution step, the powermean of all these memorized successful factors
 * is used to perturb the original random distributions, in order to generate better Cr and F
 * in the next generation.
 *
 *
 * NOTE: No memory parameter is used, since the algorithm depends on the current generation (it is recommended to use a high
 * number of generations in order to exploit this adaptive feature. Running this algorithm 1000 times with just 10 Generations 
 * consecutively on the same population will yield different results than i.e. evolving just once but for 10000 Generations!)
 * 
 * @see Sk. Minhazul Islam, S. Das, S. Ghosh, S. Roy, and P. N. Suganthan, "An Adaptive Differential Evolution Algorithm with Novel Mutation and Crossover Strategies for Global Numerical Optimization", IEEE Trans. on  Systems, Man, and Cybernetics, Part B: Cybernetics, Vol. 42, No. 2, pp.  482-500, 2012.
 *
 * @author Marcus Maertens (mmarcusx@gmail.com)
 * @author Dario Izzo (dario.izzo@gmail.com)
 */

class __PAGMO_VISIBLE mde_pbx : public base
{
public:
	mde_pbx(int = 100, double = 0.15, double = 1.5, double = 1e-30, double = 1e-30);
	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
protected:
	std::string human_readable_extra() const;
	double powermean(std::vector<double>, double) const;
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<int &>(m_gen);
		ar & const_cast<double &>(m_qperc);
		ar & const_cast<double &>(m_nexp);
		ar & const_cast<double &>(m_ftol);
		ar & const_cast<double &>(m_xtol);
		ar & m_fm;
		ar & m_crm;
		ar & m_fsuccess;
		ar & m_crsuccess;
	}
	
	// Number of generations.
	const int m_gen;

	// this vector keeps track of successfull scale factors
	mutable std::vector<double> m_fsuccess;
	
	// Current scale control factor
	mutable double m_fm;
	
	// this vector keeps track of successfull crossover probabilities
	mutable std::vector<double> m_crsuccess;

    // Current crossover control factor
	mutable double m_crm;
	
	// Control parameters
	const double m_qperc;
	const double m_nexp;
		
	// Self-adaptation strategy
	const double m_ftol;
	const double m_xtol;

};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::mde_pbx)

#endif // PAGMO_ALGORITHM_MDE_PBX_H
