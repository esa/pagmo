/*****************************************************************************
 *   Copyright (C) 2015 The PaGMO development team,                          *
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

#ifndef PAGMO_ALGORITHM_WORHP_H
#define PAGMO_ALGORITHM_WORHP_H

// WORHP uses c style boolean variables: the following definition is needed when not automatically introduced
// in worhp/C_std.h
#ifndef _Bool
#define _Bool bool
#endif

#include <worhp/worhp.h>
#include <string>

#include "worhp_cpp_wrapper/worhp_param_serialization.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"


namespace pagmo { namespace algorithm {

class __PAGMO_VISIBLE worhp: public base
{
public:
	worhp();
	void evolve(pagmo::population&) const;
	pagmo::algorithm::base_ptr clone() const;
	std::string get_name() const;

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & m_params;
	}  
	
	mutable Params m_params;
};

}} // namespaces

#endif // PAGMO_ALGORITHM_WORHP_H

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::worhp)

// We clean the unnecessary _Bool definition for future sanity
#ifdef _Bool
#undef _Bool
#endif


