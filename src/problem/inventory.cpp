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

#include "../population.h"
#include "../rng.h"
#include "../types.h"
#include "base_stochastic.h"
#include "inventory.h"

namespace pagmo { namespace problem {

inventory::inventory(int weeks,int sample_size, unsigned int seed) : base_stochastic(weeks,seed),
	m_weeks(weeks),m_sample_size(sample_size)
{
	set_lb(0.0);
	set_ub(200.0);
}

base_ptr inventory::clone() const
{
	return base_ptr(new inventory(*this));
}

void inventory::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	m_drng.seed(m_seed);
	const double c=1.0,b=1.5,h=0.1;
	double retval=0;

	for (size_t i = 0; i<m_sample_size; ++i) {
		double I=0;
		for (decision_vector::size_type j = 0; j<x.size(); ++j) {
			double d = m_drng() * 100;
			retval += c * x[j] + b * std::max<double>(d-I-x[j],0) + h * std::max<double>(I+x[j]-d,0);
			I = std::max<double>(0, I + x[j] - d);
		}
	}
	f[0] = retval / m_sample_size;
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the values vector, the weights vectors and the max weight.
 */
std::string inventory::human_readable_extra() const
{
	std::ostringstream oss;
	oss << "\n\tWeeks: " << m_weeks << '\n';
	oss << "\tSample Size: " << m_sample_size << '\n';
	oss << "\tSeed: " << m_seed << '\n';
	return oss.str();
}

std::string inventory::get_name() const
{
	return "Inventory";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::inventory)
