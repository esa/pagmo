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
#include "../src/pagmo.h"
#include "../src/util/hypervolume.h"
#include "../src/util/hv_algorithm/base.h"
#include "../src/util/hv_algorithm/hv2d.h"
#include "../src/util/hv_algorithm/hv3d.h"
#include "../src/util/hv_algorithm/hv4d.h"
#include "../src/util/hv_algorithm/wfg.h"
#include "../src/util/hv_algorithm/bf_approx.h"
#include "../src/util/hv_algorithm/bf_fpras.h"
#include "../src/util/hv_algorithm/hoy.h"
#include "../src/util/hv_algorithm/fpl.h"

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

	// Implementation of the hypervolume algorithms
	// Tested using least contributor method
	util::hv_algorithm::base_ptr bf = util::hv_algorithm::bf_approx(true, 99, 0.9, 0.99, 0.9, 0.9, 0.9, 0.9).clone();
	util::hv_algorithm::base_ptr bf_new = util::hv_algorithm::bf_approx().clone();

	// save/load the object
	{
		std::ofstream ofs("test.ar");
		boost::archive::text_oarchive oa(ofs);
		oa & bf;
	}
	{
		std::ifstream ifs("test.ar");
		boost::archive::text_iarchive ia(ifs);
		ia & bf_new;
	}

	util::hypervolume hv_7d(boost::shared_ptr<population>(new population(problem::dtlz(1, 10,7), 100)));
	util::hypervolume hv_4d(boost::shared_ptr<population>(new population(problem::dtlz(1, 10,4), 100)));
	util::hypervolume hv_3d(boost::shared_ptr<population>(new population(problem::dtlz(1, 10,3), 100)));
	util::hypervolume hv_2d(boost::shared_ptr<population>(new population(problem::dtlz(1, 10,2), 100)));
	fitness_vector nadir_7d = hv_7d.get_nadir_point(1.0);
	fitness_vector nadir_4d = hv_4d.get_nadir_point(1.0);
	fitness_vector nadir_3d = hv_3d.get_nadir_point(1.0);
	fitness_vector nadir_2d = hv_2d.get_nadir_point(1.0);

	unsigned int lc1 = hv_7d.least_contributor(nadir_7d, bf);
	unsigned int lc2 = hv_7d.least_contributor(nadir_7d, bf_new);
	if (lc1 != lc2) {
		return 1;
	}

	// Algorithms for which the compute method is available
	std::vector<util::hv_algorithm::base_ptr> compute_algs;
	std::vector<util::hv_algorithm::base_ptr> compute_algs_new;
	std::vector<std::pair<util::hypervolume, fitness_vector> > compute_input;
	compute_algs.push_back(util::hv_algorithm::bf_fpras(0.12345, 0.12345).clone());
	compute_algs_new.push_back(util::hv_algorithm::bf_fpras().clone());
	compute_input.push_back(std::make_pair(hv_7d, nadir_7d));

	compute_algs.push_back(util::hv_algorithm::wfg(4).clone());
	compute_algs_new.push_back(util::hv_algorithm::wfg().clone());
	compute_input.push_back(std::make_pair(hv_7d, nadir_7d));

	compute_algs.push_back(util::hv_algorithm::hoy().clone());
	compute_algs_new.push_back(util::hv_algorithm::hoy().clone());
	compute_input.push_back(std::make_pair(hv_3d, nadir_3d));

	compute_algs.push_back(util::hv_algorithm::hv4d().clone());
	compute_algs_new.push_back(util::hv_algorithm::hv4d().clone());
	compute_input.push_back(std::make_pair(hv_4d, nadir_4d));

	compute_algs.push_back(util::hv_algorithm::hv2d(false).clone());
	compute_algs_new.push_back(util::hv_algorithm::hv2d().clone());
	compute_input.push_back(std::make_pair(hv_2d, nadir_2d));

	compute_algs.push_back(util::hv_algorithm::hv3d(false).clone());
	compute_algs_new.push_back(util::hv_algorithm::hv3d().clone());
	compute_input.push_back(std::make_pair(hv_3d, nadir_3d));

	for(unsigned int i = 0 ; i < compute_algs.size() ; ++i) {
		// save/load the objects
		{
			std::ofstream ofs("test.ar");
			boost::archive::text_oarchive oa(ofs);
			oa & compute_algs[i];
		}
		{
			std::ifstream ifs("test.ar");
			boost::archive::text_iarchive ia(ifs);
			ia & compute_algs_new[i];
		}
		util::hypervolume &hv = compute_input[i].first;
		fitness_vector &r = compute_input[i].second;

		double hv1 = hv.compute(r, compute_algs[i]);
		double hv2 = hv.compute(r, compute_algs_new[i]);
		if (hv1 != hv2) {
			return 1;
		}
	}

	return 0;
}
