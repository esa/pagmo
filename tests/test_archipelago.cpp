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

// Test for getters and setters

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cassert>
#include "../src/pagmo.h"

using namespace pagmo;

int test_distribution_type() {
	std::vector<archipelago::distribution_type> types = {
		archipelago::distribution_type::point_to_point,
		archipelago::distribution_type::broadcast,
	};
	archipelago a;

	for(auto it = types.begin(); it != types.end(); it++) {
		a.set_distribution_type(*it);
		if(a.get_distribution_type() != *it) {
			return 1;
		}
	}

	return 0;
}

int main() {
	return test_distribution_type();
}
