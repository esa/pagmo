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

#include <iostream>
#include <iomanip>
#include "src/pagmo.h"

using namespace pagmo;

// Example in C++ of the use of PaGMO 1.1.5

int main() {
	double x_1 [16] = {10466.30234741616, 0.4134174760471763, 0.5028188179070112, 198.1961705829838, 1.548984035694917, 1.1999823933211586, 0.9999876422503512, 2.5808795876224186, -4.861037165832016, 1.0373216147703974, 0.5022120385344772, 54.908125153344876, 1.0036684306680947, 1.2570747046110966, 0.06301727892279088, 28.628205244720036};
	fitness_vector x(x_1,x_1+16);
	std::vector<kep_toolbox::planet_ptr> seq;
	seq.push_back(kep_toolbox::planet_js("callisto").clone());
	seq.push_back(kep_toolbox::planet_js("ganymede").clone());
	seq.push_back(kep_toolbox::planet_js("ganymede").clone());
	seq.push_back(kep_toolbox::planet_js("ganymede").clone());
	
	std::vector<std::vector<double> > tofs;
	std::vector<double> dumb(2);
			dumb[0] = 100;dumb[1] = 200;
			tofs.push_back(dumb);
			dumb[0] = 0.1;dumb[1] = 5;
			tofs.push_back(dumb);
			dumb[0] = 30;dumb[1] = 100;
			tofs.push_back(dumb);
			dumb[0] = 10;dumb[1] = 50;
			tofs.push_back(dumb);
	problem::mga_incipit prob(seq, kep_toolbox::epoch(7305.0), kep_toolbox::epoch(11323.0),tofs);
	fitness_vector retval = prob.objfun(x);
	std::cout << retval << std::endl;
}
