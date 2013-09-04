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

#ifndef PAGMO_PROBLEM_CON2UNCON_H
#define PAGMO_PROBLEM_CON2UNCON_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "cec2006.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// Constrained to unconstrained meta-problem
/**
 * Implements a meta-problem class that wraps some other constrained problems,
 * resulting in unconstrained problem by simply removing the constraints.
 *
 * @author Jeremie Labroquere (jeremie.labroquere@gmail.com)
 */

class __PAGMO_VISIBLE con2uncon : public base
{
public:
	//constructors
	con2uncon(const base & = cec2006(4));

	//copy constructor
	con2uncon(const con2uncon &);
	base_ptr clone() const;
	std::string get_name() const;

protected:
	std::string human_readable_extra() const;
	void objfun_impl(fitness_vector &, const decision_vector &) const;
	bool compare_fitness_impl(const fitness_vector &v_f1, const fitness_vector &v_f2) const;

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & m_original_problem;
	}
	base_ptr m_original_problem;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::con2uncon);

#endif // PAGMO_PROBLEM_con2uncon_H
