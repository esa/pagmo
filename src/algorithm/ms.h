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

#ifndef PAGMO_ALGORITHM_MS_H
#define PAGMO_ALGORITHM_MS_H

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/version.hpp>
#include <string>

#include "../config.h"
#include "../population.h"
#include "base.h"


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
	ms(const algorithm::base &, int);
	base_ptr clone() const;
	void evolve(population &) const;
	void screen_output(const bool);
	std::string get_name() const;
protected:
	std::string human_readable_extra() const;
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version){
		std::cout << "de-/serializing ms algorithm " << version << std::endl;
		ar & boost::serialization::base_object<base>(*this);
		ar & m_algorithm;
		ar & const_cast<int &>(m_starts);
		ar & m_screen_out;
	}
	boost::shared_ptr<base> m_algorithm;
	int m_starts;
	bool m_screen_out;
};

}} //namespaces

#endif // PAGMO_ALGORITHM_MS_H
