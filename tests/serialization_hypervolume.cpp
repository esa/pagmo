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
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "../src/util/hypervolume.h"

using namespace pagmo;
int main()
{
	unsigned int points_size = 10;
	unsigned int dim_size = 3;

	std::vector<fitness_vector> points(points_size, fitness_vector(dim_size, 0.0));
	for(unsigned int i = 0 ; i < points_size ; ++i) {
		for(unsigned int j = 0 ; j < dim_size ; ++j) {
			points[i][j] = i*dim_size + j;
		}
	}

	util::hypervolume_ptr hv_obj = util::hypervolume(points).clone();
	util::hypervolume_ptr hv_obj_new = util::hypervolume().clone();

	// create and open a character archive for output
	{
	std::ofstream ofs("test.ar");
	// save data to archive
	boost::archive::text_oarchive oa(ofs);
	// write class instance to archive
	oa & hv_obj;
	} // archive closes on destructor
	
	// create and open an archive for input
	{
	std::ifstream ifs("test.ar");
	// load data from the archive
	boost::archive::text_iarchive ia(ifs);
	// read class state from archive
	ia & hv_obj_new;
	} // archive closes on destructor

	for(unsigned int i = 0 ; i < points_size ; ++i) {
		for(unsigned int j = 0 ; j < dim_size ; ++j) {
			if (hv_obj->get_points()[i][j] != hv_obj_new->get_points()[i][j]) {
				std::cout << "Different points found! p_idx:" << i << " d_idx:" << j << ": " << hv_obj->get_points()[i][j] << " vs " << hv_obj_new->get_points()[i][j] << std::endl;
				return 1;
			}
		}
	}
	return 0;
}
