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

#ifndef PAGMO_ALGORITHM_BEECOLONY_H
#define PAGMO_ALGORITHM_BEECOLONY_H

#include "../config.h"
#include "base.h"


namespace pagmo { namespace algorithm {

 /*
  * The Artificial Bee Colony (ABC) algorithm is a swarm based meta-heuristic algorithm based on the behaviour of bees
  * @see http://www.scholarpedia.org/article/Artificial_bee_colony_algorithm
  * @author Andrea Mambrini (andrea.mambrini@gmail.com)
  */

class __PAGMO_VISIBLE bee_colony: public base
{
public:
	bee_colony(int gen, double onlooker_fraction=0.5, int scouts_number = 1, int limit = 20, double phi = 1.0);
	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
protected:
	std::string human_readable_extra() const;
private:
	// Number of generations.
	const int m_gen;
	const double m_onlooker_fraction;
	const int m_scouts_number;
	const int m_limit;
	const double m_phi;
};

}} //namespaces

#endif // BEECOLONY_H
