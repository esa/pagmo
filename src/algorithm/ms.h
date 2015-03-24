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

#ifndef PAGMO_ALGORITHM_MS_H
#define PAGMO_ALGORITHM_MS_H

#include <string>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"
#include "de.h"


namespace pagmo { namespace algorithm {

/// Multistart
/**
 *
 * The name of this algorithm says it all!!
 * It runs the same algorithm over and over on a random initial population. At the end, the champion will
 * keep memory of the luckiest run. The psuedo algorithm is as follow:

@verbatim
> Select a pagmo::population
> Select a pagmo::algorithm
> Store best individual
> Until termination:
> > Reset the population
> > evolve the population with the pagmo::algorithm
@endverbatim
 *
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 */

class __PAGMO_VISIBLE ms: public base
{
public:
	ms(const base & = de(), int = 1);
	ms(const ms &);
	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
	base_ptr get_algorithm() const;
	void set_algorithm(const base &);
protected:
	std::string human_readable_extra() const;
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & m_algorithm;
		ar & m_starts;
	}
	base_ptr m_algorithm;
	int m_starts;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::ms)

#endif // PAGMO_ALGORITHM_MS_H
