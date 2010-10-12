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

#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "src/pagmo.h"

using namespace pagmo;

static std::ofstream output_file("results.txt");

static const int n_proc = 48;

static std::vector<kep_toolbox::asteroid_gtoc5> gtoc_data;

struct edelbaum_comparer
{
	edelbaum_comparer(const std::vector<double> &edelbaum_dv):m_edelbaum_dv(edelbaum_dv) {}
	bool operator()(int i, int j) const
	{
		return m_edelbaum_dv[i - 1] < m_edelbaum_dv[j - 1];
	}
	const std::vector<double> &m_edelbaum_dv;
};

static inline std::vector<int> edelbaum(int i)
{
	using namespace std;
	std::vector<int> retval(7075);
	for (int j = 0; j < 7075; ++j) {
		retval[j] = j + 1;
	}
	const kep_toolbox::asteroid_gtoc5 &source_ast = gtoc_data[i - 1];
	// Calculate edelbaum.
	std::vector<double> edelbaum_dv(7075);
	kep_toolbox::array6D source_elements = source_ast.get_elements(kep_toolbox::epoch(0)), target_elements;
	for (int j = 0; j < 7075; ++j) {
		const kep_toolbox::asteroid_gtoc5 &target_ast = gtoc_data[retval[j] - 1];
		target_elements = target_ast.get_elements(kep_toolbox::epoch(0));
		const double a1 = source_elements[0], i1 = source_elements[2], W1 = source_elements[3];
		const double a2 = target_elements[0], i2 = target_elements[2], W2 = target_elements[3];
		// Circular velocity on first orbit.
		const double vc1 = std::sqrt(ASTRO_MU_SUN / a1), vc2 = std::sqrt(ASTRO_MU_SUN / a2);
		const double cos_i_rel = cos(i1)*cos(i2) + sin(i1)*sin(i2)*cos(W1)*cos(W2) + sin(i1)*sin(i2)*sin(W1)*sin(W2),
			i_rel = acos(cos_i_rel);
		edelbaum_dv[j] = sqrt(vc1 * vc1 - 2. * vc1 * vc2 * cos(boost::math::constants::pi<double>() / 2. * i_rel) + vc2 * vc2);
	}
	std::sort(retval.begin(),retval.end(),edelbaum_comparer(edelbaum_dv));
	return retval;
}

static inline std::vector<int> get_top_candidates(const std::vector<int> &visited)
{
	// Get list of best according to edelbaum.
	const std::vector<int> edelbaum_best(edelbaum(visited.back()));
	std::vector<int> retval;
	// Get top choices, excluding those already visited.
	int cur_index = 0;
	while (retval.size() < n_proc && cur_index < edelbaum_best.size()) {
		if (std::find(visited.begin(),visited.end(),edelbaum_best[cur_index]) == visited.end()) {
			retval.push_back(edelbaum_best[cur_index]);
		}
		++cur_index;
	}
	return retval;
}

static inline bool optimized_ranker(const std::pair<int,mpi_island> &a, const std::pair<int,mpi_island> &b)
{
	return (a.second.get_population().champion().x[1] + a.second.get_population().champion().x[2]) <
		(b.second.get_population().champion().x[1] + b.second.get_population().champion().x[2]);
}

static void attempt_last_fb_and_log(std::vector<double> start_mjd, std::vector<double> start_mass, std::vector<int> cur_sequence)
{
	std::cout << "Attempting final flyby\n";
	int n_segments = 5;
	algorithm::snopt algo(1000);
	const double pert_epoch = 1000, pert_nondim = 1E-1, pert_mass = 200, pert_vinf = 1000;
	std::vector<double> perturb(n_segments * 3 + 5,pert_nondim);
	perturb[0] = pert_epoch;
	perturb[1] = pert_mass;
	perturb[2] = pert_vinf;
	perturb[3] = pert_vinf;
	perturb[4] = pert_vinf;
	algorithm::mbh algo2(algo,20,perturb);
	archipelago a;
	std::cout << cur_sequence.back() << ',' << start_mass.back() << ',' << start_mjd.back() << '\n';
	for (unsigned i = 0; i < 16; ++i) {
		problem::gtoc5_self_flyby prob(n_segments,cur_sequence.back(),start_mjd.back(),start_mass.back());
		a.push_back(island(prob,algo2,1));
	}
	a.evolve(1);
	a.join();
	int i = 0;
	while (i < 16 && !a.get_island(i)->get_population().problem().feasibility_x(a.get_island(i)->get_population().champion().x)) {
		++i;
	}
	if (i < 16) {
		std::cout << "Final flyby succeeded\n";
		cur_sequence.push_back(cur_sequence.back());
		start_mass.push_back(a.get_island(i)->get_population().champion().x[1]);
		start_mjd.push_back(a.get_island(i)->get_population().champion().x[0] + start_mjd.back());
	} else {
		std::cout << "Final flyby failed\n";
	}
	output_file << cur_sequence << '\n';
	output_file << start_mass << '\n';
	output_file << start_mjd << '\n';
}

static void recursive_step(std::vector<double> start_mjd, std::vector<double> start_mass, std::vector<int> cur_sequence)
{
	std::cout << "Current sequence: " << cur_sequence << '\n';

	//1 - Locate next plausible targets (max n_proc)
	std::vector<int> top_choices = get_top_candidates(cur_sequence);
	if (top_choices.size() == 0) {
		attempt_last_fb_and_log(start_mjd,start_mass,cur_sequence);
		return;
	}

	std::cout << "Top choices Edelbaum: " << top_choices << '\n';

	//2 - Optimize all trajectories to plausible targets
	archipelago a;
	// create the algorithms
	int n_segments = 5;
	algorithm::snopt algo(1000);
	const double pert_epoch = 1000, pert_nondim = 1E-1, pert_mass = 200, pert_vinf = 1000;
	std::vector<double> perturb(n_segments * 6 + 8,pert_nondim);
	perturb[0] = pert_epoch;
	perturb[1] = pert_epoch;
	perturb[2] = pert_epoch;
	perturb[3] = pert_mass;
	perturb[4] = pert_mass;
	perturb[5] = pert_vinf;
	perturb[6] = pert_vinf;
	perturb[7] = pert_vinf;
	algorithm::mbh algo2(algo,20,perturb);
	algorithm::ms algo3(algo2,5);
	for (unsigned i = 0; i < top_choices.size(); ++i) {
		problem::gtoc5_flyby prob(n_segments,cur_sequence.back(),cur_sequence.back(),top_choices[i],start_mjd.back(),start_mass.back(),problem::gtoc5_flyby::MASS);
		a.push_back(mpi_island(prob,algo3,5));
	}
// 	a.evolve(1);
// 	a.join();
	std::vector<std::pair<int,mpi_island> > bah;
	for (unsigned i = 0; i < top_choices.size(); ++i) {
		// Strip unfeasbile solutions
		if (a.get_island(i)->get_population().problem().feasibility_x(a.get_island(i)->get_population().champion().x)) {
			bah.push_back(std::make_pair(top_choices[i],*dynamic_cast<mpi_island *>(a.get_island(i).get())));
		}
	}
	// If there are no feasbile solutions DEAD END!
	if (bah.size() == 0) {
		attempt_last_fb_and_log(start_mjd,start_mass,cur_sequence);
		return;
	}

	//3 - Rank the top choices according to the optimized results
	std::sort(bah.begin(),bah.end(),optimized_ranker);

	std::cout << "Top choices after optimization: [";
	for (unsigned i = 0; i < bah.size(); ++i) {
		if (i != bah.size() - 1) {
			std::cout << bah[i].first << ',';
		}
	}
	std::cout << "]\n";

	// 4 - For each feasible solution, check DV/T limits and either start new recursive_step or stop
	for (unsigned i = 0; i < bah.size(); ++i) {
		const double final_mass = bah[i].second.get_population().champion().x[4];
		const double final_time = bah[i].second.get_population().champion().x[0] + bah[i].second.get_population().champion().x[1] + bah[i].second.get_population().champion().x[2] - start_mjd[0];
		if (final_mass < 350 || final_time > 365.25 * 16) {
			attempt_last_fb_and_log(start_mjd,start_mass,cur_sequence);
			return;
		} else {
			start_mjd.push_back(final_time + start_mjd[0]);
			start_mass.push_back(final_mass);
			cur_sequence.push_back(bah[i].first);
			recursive_step(start_mjd,start_mass,cur_sequence);
		}
	}
}

int main()
{
	mpi_environment env;

	// Init asteroid data.
	for (int i = 1; i <= 7075; ++i) {
		gtoc_data.push_back(kep_toolbox::asteroid_gtoc5(i));
	}

	// Starting asteroid id.
	std::vector<int> start_sequence(1,5249);
	// Starting mjd.
	std::vector<double> start_mjd(1,58890.282931514834 + 356.5346936499468);
	std::vector<double> start_mass(1,3985.4461077866031 - 40);

	recursive_step(start_mjd,start_mass,start_sequence);
	
// 	int n_segments = 10;
// 
// 	algorithm::snopt algo(1000);
// 
// 	const double pert_epoch = 1000;
// 	const double pert_nondim = 1E-1;
// 	const double pert_mass = 200;
// 	const double pert_vinf = 1000;
// 	std::vector<double> perturb(n_segments * 3 + 6,pert_nondim);
// 	perturb[0] = pert_epoch;
// 	perturb[1] = pert_epoch;
// 	perturb[2] = pert_mass;
// 	perturb[3] = pert_vinf;
// 	perturb[4] = pert_vinf;
// 	perturb[5] = pert_vinf;
// 
// 	algorithm::mbh algo2(algo,20,perturb);
// 	algorithm::ms algo3(algo2,5);
// 	//algo2.screen_output(true);
// 
// 	// Outer loop.
// 	int cur_id = 1;
// 	for (int i = 1; i <= 25; ++i) {
// 		archipelago a;
// 		for (int j = 0; j < 283; ++j) {
// 			a.push_back(mpi_island(problem::gtoc5_launch(n_segments,cur_id + j,problem::gtoc5_launch::TIME),algo3,5));
// 		}
// 		a.evolve(1);
// 		a.join();
// 		std::ofstream ofs((std::string("time_") + boost::lexical_cast<std::string>(cur_id)).c_str());
// 		boost::archive::text_oarchive oa(ofs);
// 		oa << static_cast<const archipelago &>(a);
// 		cur_id += 283;
// 	}
// 
// 	cur_id = 1;
// 	for (int i = 1; i <= 25; ++i) {
// 		archipelago a;
// 		for (int j = 0; j < 283; ++j) {
// 			a.push_back(mpi_island(problem::gtoc5_launch(n_segments,cur_id + j,problem::gtoc5_launch::MASS),algo3,5));
// 		}
// 		a.evolve(1);
// 		a.join();
// 		std::ofstream ofs((std::string("mass_") + boost::lexical_cast<std::string>(cur_id)).c_str());
// 		boost::archive::text_oarchive oa(ofs);
// 		oa << static_cast<const archipelago &>(a);
// 		cur_id += 283;
// 	}
}
