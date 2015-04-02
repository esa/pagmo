/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 3 of the License, or       *
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
#include <string>
#include <keplerian_toolbox/keplerian_toolbox.h>

#include "../src/algorithm/snopt.h"
#include "../src/island.h"
#include "../src/problem/gtoc_2.h"

using namespace pagmo;
using namespace kep_toolbox;

/*
In this example we instantiate a problenm of the type gtoc_2 using the asteroid sequence which was used by the Turin team
to win the 2nd edition of the global trajectory optimization competition (http://www.esa.int/gsp/ACT/mad/op/GTOC/index.htm), later we define
a starting guess containing the epochs taken from the winning solution and we solve the problem.
*/

int main()
{
	//Instantiate the problem with default 10 segments
	problem::gtoc_2 prob(815,300,110,47);

	//Create a population containing a single random individual
	population pop(prob,1);

	//Set the Turin solution epochs
 	decision_vector tmp = pop.get_individual(0).cur_x;
	tmp[0] = 59870; tmp[1] = 60283 - 59870; tmp[2] = 60373 - 60283;
 	tmp[3] = 61979 - 60373; tmp[4] = 62069 - 61979; tmp[5] = 62647 - 62069;
	tmp[6] = 62737 - 62647; tmp [7] = 63196 - 62737;
	tmp[8] = 1400; tmp[9] = 1200; tmp[10]= 1100; tmp[11] = 1000;
 	pop.set_x(0, tmp);

	//Instantiate the algorithm
	algorithm::snopt algo(1000,1E-9,1E-9);
	algo.set_screen_output(true);
	
	//Create the island  
	island isl(algo,pop);

	//Solve the problem
	isl.evolve();

	std::cout << isl << '\n';

	return 0;
}
