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

#include <worhp/worhp.h>
#include <string>
#include <map>

#include "../serialization.h"
#include "worhp_cpp_wrapper/worhp_param_serialization.h"
#include "base.h"

namespace pagmo {

// Forward declaration
class population;

namespace algorithm {

class __PAGMO_VISIBLE worhp: public base
{
public:
	worhp(const int iter=100, const double feas=1e-10, const double opt=1e-4, const bool screen_output=true);
	void evolve(pagmo::population&) const;
	pagmo::algorithm::base_ptr clone() const;
	std::string get_name() const;

	void set_param(const std::string, const double);
	double get_param(const std::string name) const;
	std::vector<std::string> get_available_parameters() const;

protected:
	std::string human_readable_extra() const;

private:
	void define_param_map() ;

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & m_params;
	}
	// Structure containing all the WORHP parameters, see WORHP documentation
	Params m_params;
	// A map between the WORHP parameters available to the user for geting and setting and integres 1...N
	std::map<std::string, int> m_param_map;
};

}} // namespaces

#endif // PAGMO_ALGORITHM_WORHP_H

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::worhp)
