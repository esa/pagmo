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

#include <boost/numeric/conversion/cast.hpp>
#include <boost/random/uniform_int.hpp>
#include <climits>
#include <cmath>
#include <cstddef>
#include <exception>
#include <sstream>
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../types.h"
#include "base.h"
#include "ihs.h"

namespace pagmo
{
namespace algorithm {

/// Simple constructor.
/**
 * Initialises the parameters of the algorithms to "standard" values: phmcr = 0.85, ppar_min = 0.35, ppar_max = 0.99,
 * bw_min = 1E-5 and bw_max = 1.
 *
 * @param[in] gen number of generations.
 * @throws value_error if gen is not positive.
 */
ihs::ihs(int gen):base(),m_gen(gen),m_phmcr(0.85),m_ppar_min(0.35),m_ppar_max(0.99),m_bw_min(1E-5),m_bw_max(1)
{
	boost::numeric_cast<std::size_t>(gen);
	if (!m_gen) {
		pagmo_throw(value_error,"number of generations must be positive");
	}
}

/// Advanced constructor.
/**
 * Allows to specify in detail the parameters of the algorithm.
 *
 * @param[in] gen number of generations.
 * @param[in] phmcr rate of choosing from memory.
 * @param[in] ppar_min minimum pitch adjustment rate.
 * @param[in] ppar_max maximum pitch adjustment rate.
 * @param[in] bw_min minimum distance bandwidth.
 * @param[in] bw_max maximum distance bandwidth.
 * @throws value_error if gen is not positive, phmcr is not in the ]0,1[ interval, ppar min/max are not in the ]0,1[ interval,
 * min/max quantities are less than/greater than max/min quantities.
 */
ihs::ihs(int gen, const double &phmcr, const double &ppar_min, const double &ppar_max, const double &bw_min, const double &bw_max):
	base(),m_gen(gen),m_phmcr(phmcr),m_ppar_min(ppar_min),m_ppar_max(ppar_max),m_bw_min(bw_min),m_bw_max(bw_max)
{
	boost::numeric_cast<std::size_t>(gen);
	if (!m_gen) {
		pagmo_throw(value_error,"number of generations must be positive");
	}
	if (phmcr > 1 || phmcr < 0 || ppar_min > 1 || ppar_min < 0 || ppar_max > 1 || ppar_max < 0) {
		pagmo_throw(value_error,"probability of choosing from memory and pitch adjustment rates must be in the [0,1] range");
	}
	if (ppar_min > ppar_max) {
		pagmo_throw(value_error,"minimum pitch adjustment rate must not be greater than maximum pitch adjustment rate");
	}
	if (bw_min <= 0 || bw_max < bw_min) {
		pagmo_throw(value_error,"bandwidth values must be positive, and minimum bandwidth must not be greater than maximum bandwidth");
	}
}

/// Clone method.
base_ptr ihs::clone() const
{
	return base_ptr(new ihs(*this));
}

/// Evolve implementation.
/**
 * Run the IHS algorithm for the number of generations specified in the constructors. Within each call of this method,
 * the ppar and bw parameters will be varied between maximum and minimum values according to the IHS schedule.
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void ihs::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type prob_dimension = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type pop_size = pop.size();
	// Go out if there is nothing to do.
	if (pop_size == 0 || m_gen == 0) {
		return;
	}
	// Init temporary decision vector, and lower-upper bounds difference vector.
	m_tmp_x.resize(prob_dimension);
	m_lu_diff.resize(prob_dimension);
	for (problem::base::size_type i = 0; i < prob_dimension; ++i) {
		m_lu_diff[i] = ub[i] - lb[i];
	}
	// Temporary fitness vector.
	m_tmp_f.resize(prob.get_f_dimension());
	const double c = std::log(m_bw_min/m_bw_max) / m_gen;
	for (std::size_t g = 0; g < m_gen; ++g) {
		const double ppar_cur = m_ppar_min + ((m_ppar_max - m_ppar_min) * g) / m_gen, bw_cur = m_bw_max * std::exp(c * g);
		for (problem::base::size_type i = 0; i < prob_dimension; ++i) {
			const double next_rn = m_drng();
			if (m_drng() < m_phmcr) {
				// With random probability, tmp's i-th chromosome element is the one from a randomly chosen individual.
				m_tmp_x[i] = pop.get_individual((population::size_type)(next_rn * pop_size)).get<0>()[i];
				if (m_drng() < ppar_cur) {
					if (i < prob_dimension - prob_i_dimension) {
						// Handle the continuous part of the problem.
						// Randomly, add or subtract pitch from the current chromosome element.
						const double next_next_rn = m_drng();
						if (m_drng() > .5) {
							m_tmp_x[i] += next_next_rn * bw_cur * m_lu_diff[i];
						} else {
							m_tmp_x[i] -= next_next_rn * bw_cur * m_lu_diff[i];
						}
					} else {
						// Integer part of the problem.
						if (m_tmp_x[i] == INT_MIN || m_tmp_x[i] == INT_MAX) {
							pagmo_throw(std::overflow_error,"integer overflow in harmony search algortihm");
						}
						if (m_drng() > .5) {
							m_tmp_x[i] += 1;
						} else {
							m_tmp_x[i] -= 1;
						}
					}
					// Handle the case in which we addded or subtracted too much and ended up out
					// of boundaries.
					if (m_tmp_x[i] > ub[i]) {
						m_tmp_x[i] = ub[i];
					} else if (m_tmp_x[i] < lb[i]) {
						m_tmp_x[i] = lb[i];
					}
				}
			} else {
				// Pick randomly within the bounds.
				if (i < prob_dimension - prob_i_dimension) {
					// Continuous.
					m_tmp_x[i] = lb[i] + next_rn * m_lu_diff[i];
				} else {
					// Integral.
					m_tmp_x[i] = boost::uniform_int<int>(lb[i],ub[i])(m_urng);
				}
			}
		}
		const population::size_type worst_idx = pop.get_worst_idx();
		prob.objfun(m_tmp_f,m_tmp_x);
		if (prob.compare_f(m_tmp_f,pop.get_individual(worst_idx).get<2>())) {
			pop.set_x(worst_idx,m_tmp_x);
		}
	}
}

/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string ihs::human_readable_extra() const
{
	std::ostringstream s;
	s << "\tgenerations:\t" << m_gen << '\n';
	s << "\tphmcr:\t\t" << m_phmcr << '\n';
	s << "\tppar_min:\t" << m_ppar_min << '\n';
	s << "\tppar_max:\t" << m_ppar_max << '\n';
	s << "\tbw_min:\t\t" << m_bw_min << '\n';
	s << "\tbw_max:\t\t" << m_bw_max << '\n';
	return s.str();
}

}
}
