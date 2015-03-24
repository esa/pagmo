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

#ifndef PAGMO_ALGORITHM_IHS_H
#define PAGMO_ALGORITHM_IHS_H

#include <cstddef>
#include <iostream>
#include <string>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"

namespace pagmo { namespace algorithm {

// TODO: automatic reshuffling when the population becomes too homogeneous.
/// Improved harmony search algorithm.
/**
 * Harmony search (HS) is a metaheuristic algorithm mimicking the improvisation process of musicians. In the process, each musician (i.e., each variable)
 * plays (i.e., generates) a note (i.e., a value) for finding a best harmony (i.e., the global optimum) all together.
 *
 * This code implements the so-called improved harmony search algorithm (IHS), in which the probability of picking the variables from
 * the decision vector and the amount of mutation to which they are subject vary respectively linearly and exponentially within each call
 * of the evolve() method.
 *
 * In this algorithm the number of objective function evaluations is equal to the number of generations. All the individuals in the input population participate
 * in the evolution. A new individual is generated at every iteration, substituting the current worst individual of the population if better. This algorithm
 * will use the comparison methods provided by the problem in order to rank individuals.
 *
 * This algorithm is suitable for continuous, constrained, mixed-integer and multi-objective optimisation.
 *
 * @see http://en.wikipedia.org/wiki/Harmony_search for an introduction on harmony search.
 * @see http://dx.doi.org/10.1016/j.amc.2006.11.033 for the paper that introduces and explains improved harmony search.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 * @author Dario Izzo (dario.izzo@googlemail.com)
 */
class __PAGMO_VISIBLE ihs: public base
{
	public:
		ihs(int gen = 1, const double &phmcr = 0.85, const double &ppar_min = 0.35, const double &ppar_max = 0.99,
			const double &bw_min = 1E-5, const double &bw_max = 1);
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
			ar & const_cast<double &>(m_phmcr);
			ar & const_cast<double &>(m_ppar_min);
			ar & const_cast<double &>(m_ppar_max);
			ar & const_cast<double &>(m_bw_min);
			ar & const_cast<double &>(m_bw_max);
		}
		// Number of generations.
		const std::size_t		m_gen;
		// Rate of choosing from memory (i.e., from population).
		const double			m_phmcr;
		// Minimum pitch adjustment rate.
		const double			m_ppar_min;
		// Maximum pitch adjustment rate.
		const double			m_ppar_max;
		// Mininum distance bandwidth.
		const double			m_bw_min;
		// Maximum distance bandwidth.
		const double			m_bw_max;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::ihs)

#endif
