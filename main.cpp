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

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "src/pagmo.h"
#include "src/keplerian_toolbox/keplerian_toolbox.h"

using namespace pagmo;
using namespace kep_toolbox;

int main()
{

	problem::gtoc_2 prob(815,300,110,47,10);

	population pop(prob,1);
	decision_vector tmp = pop.get_individual(0).cur_x;
	tmp[0] = 59870; tmp[1] = 60283 - 59870;
	tmp[3] = 61979 - 60373; tmp[5] = 62647 - 62069;
	tmp [7] = 63196 - 62737;
	tmp[8] = 1300; tmp[9] = 1100; tmp[10]= 900; tmp[11] = 700;
	pop.set_x(0, tmp);
	algorithm::snopt algo(1000,1E-9,1E-9);
	algo.screen_output(true);
	
	island isl(pop,algo);

	//std::cout << prob << std::endl;

	isl.evolve();isl.join();

	std::cout << prob.pretty(isl.get_population().champion().x) << std::endl;

	return 0;
}
