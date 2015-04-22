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
#include <string>
#include <vector>

#include "../worhp.h"

namespace pagmo { namespace algorithm {

void worhp::define_param_map() {
	m_param_map["AcceptTolFeas"] = 1;
	m_param_map["AcceptTolOpti"] = 2;
	m_param_map["TolFeas"] = 3;
	m_param_map["TolComp"] = 4;
	m_param_map["TolOpti"] = 5;
	m_param_map["NLPprint"] = 6;
	m_param_map["InitialLMest"] = 7;
	m_param_map["MaxIter"] = 8;
}

void worhp::set_param(const std::string name, const double value)
{
	// We check that the key exists
	auto el = m_param_map.find(name);
	if (el == m_param_map.end()) {
	    pagmo_throw(value_error,"Unknown parameter name, cannot set it.");
	}
	// According to the input we set the appropriate parameter
	switch ( m_param_map[name] ) {
		case 1:
			m_params.AcceptTolFeas = value;
			break;
		case 2:
			m_params.AcceptTolOpti = value;
			break;
		case 3:
			m_params.TolFeas = value;
			break;
		case 4:
			m_params.TolComp = value;
			break;
		case 5:
			m_params.TolOpti = value;
			break;
		case 6:
			m_params.NLPprint = value;
			break;
		case 7:
			m_params.InitialLMest = value;
			break;
		case 8:
			m_params.MaxIter = value;
			break;
		default:
			pagmo_throw(value_error,"You should not be here!!");
	}
}

double worhp::get_param(const std::string name) const
{
	// We check that the key exists
	auto el = m_param_map.find(name);
	if (el == m_param_map.end()) {
		pagmo_throw(value_error,"Unknown parameter name, cannot get it.");
	}

	// We make a copy for const correctness
	std::map<std::string, int> map_copy = m_param_map;

	double retval = 0;

	// According to the input we return the appropriate parameter
	switch ( map_copy[name] ) {
	case 1:
		retval = m_params.AcceptTolFeas;
		break;
	case 2:
		retval = m_params.AcceptTolOpti;
		break;
	case 3:
		retval = m_params.TolFeas;
		break;
	case 4:
		retval = m_params.TolComp;
		break;
	case 5:
		retval = m_params.TolOpti;
		break;
	case 6:
		retval = m_params.NLPprint;
		break;
	case 7:
		retval = m_params.InitialLMest;
		break;
	case 8:
		retval = m_params.MaxIter;
		break;
	default:
		pagmo_throw(value_error,"You should not be here!!");
	}
	return retval;
}

std::vector<std::string> worhp::get_available_parameters() const
{
	std::vector<std::string> retval;
	for(auto it = m_param_map.begin(); it != m_param_map.end(); ++it) {
		retval.push_back(it->first);
	}
	return retval;
}

}} // namespaces
