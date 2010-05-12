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

#include "../population.h"
#include "../rng.h"
#include "../types.h"
#include "base.h"
#include "inventory.h"

namespace pagmo
{
namespace problem {

inventory::inventory(int sample_size):base(1),m_seed(rng_generator::get<rng_uint32>()()),
	m_sample_size(sample_size),m_drng(m_seed)
{
	set_lb(0.0);
	set_ub(100.0);
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
		double d = m_drng() * 100;
		retval += c * x[0] + b * std::max<double>(d-x[0],0.0) + h * std::max<double>(x[0] - d, 0.0);
	}
	f[0] = retval / m_sample_size;
}

bool inventory::equality_operator_extra(const base &other) const
{
	return false;
	/*return (m_d == dynamic_cast<const inventory *>(&other)->m_d && m_sample_size ==
		dynamic_cast<const inventory *>(&other)->m_sample_size);*/
}

void inventory::pre_evolution(population &pop) const
{
	m_seed = rng_generator::get<rng_uint32>()();
// 	for (size_t i = 0; i<m_sample_size; ++i) {
// 		m_d[i] = m_drng() * 100;
// 	}
	//Re-evaluate the population with respect to the new seed (Internal Sampling Method)
	for (population::size_type i = 0; i<pop.size(); ++i) {
		pop.set_x(i,pop.get_individual(i).cur_x);
	}
}

}
}
