/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
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

#ifndef PAGMO_PROBLEM_DECOMPOSITION_H
#define PAGMO_PROBLEM_DECOMPOSITION_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "base.h"
#include "ackley.h"

namespace pagmo{ namespace problem {

/// Decomposition meta-problem
/**
 * Implements a meta-problem class resulting in a decomposed version
 * of the multi-objective input problem, i.e. a single-objective problem having as fitness function 
 * a convex combination of the original fitness functions.
 *
 * @author Andrea Mambrini (andrea.mambrni@gmail.com)
 */

class __PAGMO_VISIBLE decomposition : public base
{
	public:
		//constructor
		decomposition(const base & = ackley(1), const std::vector<double> & = std::vector<double>(1,1));  //TODO: fix default vect with the proper size
		
		//copy constructor
		decomposition(const decomposition &);
		base_ptr clone() const;
		std::string get_name() const;
		
	protected:
		std::string human_readable_extra() const;
		void objfun_impl(fitness_vector &, const decision_vector &) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & m_original_problem;
			ar & m_weight;
		}
		base_ptr m_original_problem;
		std::vector<double> m_weight;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::decomposition);

#endif // PAGMO_PROBLEM_DECOMPOSITION_H
