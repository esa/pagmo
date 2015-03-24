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

#include <boost/numeric/conversion/cast.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <climits>
#include <cmath>
#include <cstddef>
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

/// Constructor.
/**
 * Allows to specify in detail the parameters of the algorithm.
 *
 * @param[in] iterations number of iterations.
 * @param[in] phmcr rate of choosing from memory.
 * @param[in] ppar_min minimum pitch adjustment rate.
 * @param[in] ppar_max maximum pitch adjustment rate.
 * @param[in] bw_min minimum distance bandwidth.
 * @param[in] bw_max maximum distance bandwidth.
 * @throws value_error if phmcr is not in the ]0,1[ interval, ppar min/max are not in the ]0,1[ interval,
 * min/max quantities are less than/greater than max/min quantities.
 */
ihs::ihs(int iterations, const double &phmcr, const double &ppar_min, const double &ppar_max, const double &bw_min, const double &bw_max):
	base(),m_gen(boost::numeric_cast<std::size_t>(iterations)),m_phmcr(phmcr),m_ppar_min(ppar_min),m_ppar_max(ppar_max),m_bw_min(bw_min),m_bw_max(bw_max)
{
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

base_ptr ihs::clone() const
{
	return base_ptr(new ihs(*this));
}

void ihs::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type prob_dimension = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type pop_size = pop.size();
	// Get out if there is nothing to do.
	if (pop_size == 0 || m_gen == 0) {
		return;
	}
	decision_vector lu_diff(prob_dimension);
	for (problem::base::size_type i = 0; i < prob_dimension; ++i) {
		lu_diff[i] = ub[i] - lb[i];
	}
	// Int distribution to be used when picking random individuals.
	boost::uniform_int<population::size_type> uni_int(0,pop_size - 1);
	const double c = std::log(m_bw_min/m_bw_max) / m_gen;
	// Temporary individual used during evolution.
	population::individual_type tmp;
	tmp.cur_x.resize(prob_dimension);
	tmp.cur_f.resize(prob.get_f_dimension());
	tmp.cur_c.resize(prob.get_c_dimension());
	for (std::size_t g = 0; g < m_gen; ++g) {
		const double ppar_cur = m_ppar_min + ((m_ppar_max - m_ppar_min) * g) / m_gen, bw_cur = m_bw_max * std::exp(c * g);
		// Continuous part.
		for (problem::base::size_type i = 0; i < prob_dimension - prob_i_dimension; ++i) {
			if (m_drng() < m_phmcr) {
				// tmp's i-th chromosome element is the one from a randomly chosen individual.
				tmp.cur_x[i] = pop.get_individual(uni_int(m_urng)).cur_x[i];
				// Do pitch adjustment with ppar_cur probability.
				if (m_drng() < ppar_cur) {
					// Randomly, add or subtract pitch from the current chromosome element.
					if (m_drng() > .5) {
						tmp.cur_x[i] += m_drng() * bw_cur * lu_diff[i];
					} else {
						tmp.cur_x[i] -= m_drng() * bw_cur * lu_diff[i];
					}
					// Handle the case in which we added or subtracted too much and ended up out
					// of boundaries.
					if (tmp.cur_x[i] > ub[i]) {
						tmp.cur_x[i] = boost::uniform_real<double>(lb[i],ub[i])(m_drng);
					} else if (tmp.cur_x[i] < lb[i]) {
						tmp.cur_x[i] = boost::uniform_real<double>(lb[i],ub[i])(m_drng);
					}
				}
			} else {
				// Pick randomly within the bounds.
				tmp.cur_x[i] = boost::uniform_real<double>(lb[i],ub[i])(m_drng);
			}
		}

		//Integer Part
		for (problem::base::size_type i = prob_dimension - prob_i_dimension; i < prob_dimension; ++i) {
			if (m_drng() < m_phmcr) {
				tmp.cur_x[i] = pop.get_individual(uni_int(m_urng)).cur_x[i];
				if (m_drng() < ppar_cur) {
					if (m_drng() > .5) {
						tmp.cur_x[i] += double_to_int::convert(m_drng() * bw_cur * lu_diff[i]);
					} else {
						tmp.cur_x[i] -= double_to_int::convert(m_drng() * bw_cur * lu_diff[i]);
					}
					// Wrap over in case we went past the bounds.
					if (tmp.cur_x[i] > ub[i]) {
						tmp.cur_x[i] = lb[i] + double_to_int::convert(tmp.cur_x[i] - ub[i]) % static_cast<int>(lu_diff[i]);
					} else if (tmp.cur_x[i] < lb[i]) {
						tmp.cur_x[i] = ub[i] - double_to_int::convert(lb[i] - tmp.cur_x[i]) % static_cast<int>(lu_diff[i]);
					}
				}
			} else {
				// Pick randomly within the bounds.
				tmp.cur_x[i] = boost::uniform_int<int>(lb[i],ub[i])(m_urng);
			}
		}
		// And we push him back
		pop.push_back(tmp.cur_x);
		// We locate the worst individual.
		const population::size_type worst_idx = pop.get_worst_idx();
		// And we get rid of him :)
		pop.erase(worst_idx);
	}
}


/// Algorithm name
std::string ihs::get_name() const
{
	return "Improved Harmony Search";
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string ihs::human_readable_extra() const
{
	std::ostringstream s;
	s << "iter:" << m_gen << ' ';
	s << "phmcr:" << m_phmcr << ' ';
	s << "ppar_min:" << m_ppar_min << ' ';
	s << "ppar_max:" << m_ppar_max << ' ';
	s << "bw_min:" << m_bw_min << ' ';
	s << "bw_max:" << m_bw_max << ' ';
	return s.str();
}

}
}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::ihs)
