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

static std::ofstream output_file("results_true.txt");

static const int n_proc = 16;

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

static inline bool optimized_ranker(const std::pair<int,island> &a, const std::pair<int,island> &b)
{
	return (a.second.get_population().champion().x[1] + a.second.get_population().champion().x[2]) <
		(b.second.get_population().champion().x[1] + b.second.get_population().champion().x[2]);
}

static void attempt_last_fb_and_log(std::vector<double> mjd_vector, std::vector<double> mass_vector, std::vector<int> seq_vector)
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
	algorithm::mbh algo2(algo,10,perturb);
	archipelago a;
	std::cout << seq_vector.back() << ',' << mass_vector.back() << ',' << mjd_vector.back() << '\n';
	for (unsigned i = 0; i < 16; ++i) {
		problem::gtoc5_self_flyby prob(n_segments,seq_vector.back(),mjd_vector.back(),mass_vector.back());
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
		seq_vector.push_back(seq_vector.back());
		mass_vector.push_back(a.get_island(i)->get_population().champion().x[1]);
		mjd_vector.push_back(a.get_island(i)->get_population().champion().x[0] + mjd_vector.back());
	} else {
		std::cout << "Final flyby failed\n";
	}
	output_file << seq_vector << '\n';
	output_file << mass_vector << '\n';
	output_file << mjd_vector << '\n';
	output_file.flush();
}

static void recursive_step(std::vector<double> mjd_vector, std::vector<double> mass_vector, std::vector<int> seq_vector)
{
	std::cout << "Current sequence: " << seq_vector << '\n';
	std::cout << "Current masses: " << mass_vector << '\n';
	std::cout << "Current epochs: [0,";
	for (int i=1; i< mjd_vector.size(); ++i){
		std::cout << (mjd_vector[i] - mjd_vector[i-1]) / 365.25 << ',';
	}
	std::cout << "]\n";

	//1 - Locate next n_proc plausible targets according to Edelbaum criteria
	std::vector<int> top_choices = get_top_candidates(seq_vector);
	std::cout << "Top choices Edelbaum: " << top_choices << '\n';
	if (top_choices.size() == 0) {
		attempt_last_fb_and_log(mjd_vector,mass_vector,seq_vector);
		return;
	}

	//2 - Optimize globally all trajectories to plausible targets
	archipelago a;
	// create the algorithms
	int n_segments = 5;
	algorithm::snopt algo(1000);
	
	//algo.screen_output(true);
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
		problem::gtoc5_flyby prob(n_segments,seq_vector.back(),seq_vector.back(),top_choices[i],mjd_vector.back(),mass_vector.back(),problem::gtoc5_flyby::TIME);
		a.push_back(island(prob,algo3,5));
	}
 	a.evolve(1);
 	a.join();
	std::vector<std::pair<int,island> > bah;
	for (unsigned i = 0; i < top_choices.size(); ++i) {
		// Strip unfeasbile solutions
		if (a.get_island(i)->get_population().problem().feasibility_x(a.get_island(i)->get_population().champion().x)) {
			bah.push_back(std::make_pair(top_choices[i],*dynamic_cast<island *>(a.get_island(i).get())));
		}
	}
	// If there are no feasible solutions DEAD END!
	if (bah.size() == 0) {
		attempt_last_fb_and_log(mjd_vector,mass_vector,seq_vector);
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
		const double final_mass = bah[i].second.get_population().champion().x[4] - 40;
		const double final_mjd = bah[i].second.get_population().champion().x[0] + bah[i].second.get_population().champion().x[1] + bah[i].second.get_population().champion().x[2];
		if (final_mass < 350 || (final_mjd - mjd_vector[0]) > 365.25 * 16) {
			attempt_last_fb_and_log(mjd_vector,mass_vector,seq_vector);
			return;
		} else {
			std::vector<double> mjd_vector_copy(mjd_vector), mass_vector_copy(mass_vector);
			std::vector<int> seq_vector_copy(seq_vector);
			// Push back the next leg's data.
			mjd_vector_copy.push_back(final_mjd);
			mass_vector_copy.push_back(final_mass);
			seq_vector_copy.push_back(bah[i].first);
			recursive_step(mjd_vector_copy,mass_vector_copy,seq_vector_copy);
		}
	}
}

int main()
{
	//mpi_environment env;

	// Init asteroid data.
	for (int i = 1; i <= 7075; ++i) {
		gtoc_data.push_back(kep_toolbox::asteroid_gtoc5(i));
	}
	// Earth
	std::vector<int> start_sequence(1,7076);
	std::vector<double> mjd_vector(1,58890.282931514834);
	std::vector<double> mass_vector(1,4000);

	// First Asteroid
	start_sequence.push_back(5249);
	mjd_vector.push_back(58890.282931514834 + 356.5346936499468);
	mass_vector.push_back(3985.4461077866031 - 40);

	//Start Recursion
	recursive_step(mjd_vector,mass_vector,start_sequence);
}
