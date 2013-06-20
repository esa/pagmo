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
#include <fstream>
#include <sstream>
#include <string>
#include "../src/pagmo.h"
#include "../src/util/hypervolumes.h"

using namespace pagmo;

int run_hypervolume_tests(std::ifstream &input, std::ofstream &output, std::string &method_name, double eps)
{
	int T, dim, s;
	double ans;
	input >> T;
	int OK_counter=0;
	for(int t=0 ; t < T ; ++t) {
		input >> dim >> s;
		std::vector<fitness_vector> ps(s, fitness_vector(dim, 0.0));
		fitness_vector r(dim, 0.0);
		for(int d = 0 ; d < dim ; ++d) {
			input >> r[d];
		}
		for(int i = 0 ; i < s ; ++i) {
			for(int d = 0 ; d < dim ; ++d) {
				input >> ps[i][d];
			}
		}
		input >> ans;

		double hypvol;
		if (method_name == "lebmeasure") {
			hypvol = pagmo::util::hypervolumes::lebmeasure(ps, r);
		} else if (method_name == "optimal2d") {
			hypvol = pagmo::util::hypervolumes::optimal2d(ps, r);
		} else {
			output << "Unknown method (" << method_name << ") .. exiting";
			exit(1);
		}

		if (fabs(hypvol-ans) < eps)
		{
			++OK_counter;
		} else {
			output << "\n Error in test " << t << ": " << hypvol << " != " << ans;
		}
	}
	output << std::endl<< " " << OK_counter << "/" << T << " passed\n";
	return (OK_counter < T ? 1 : 0);
}

int main()
{
std::string line;
std::string input_data_dir("hypervolume_testing_data/");
std::string input_data_testcases_dir(input_data_dir + "testcases/");
std::ifstream ifs((input_data_dir + "testcases_list").c_str());
std::ofstream output("hypervolume_testing_report");
output.precision(15);

int test_result = 0;
if (ifs.is_open()) {
	while(ifs.good()) {
		std::string line;
		getline(ifs,line);
		if (line == "" || line[0] == '#')
			continue;
		std::stringstream ss (line);
		std::string method_name;
		std::string test_name;
		double eps;
		ss >> method_name;
		ss >> test_name;
		ss >> eps;
		std::ifstream input((input_data_testcases_dir + test_name).c_str());
		if (input.is_open()){
			output << "Test " << test_name << " (eps:" << eps << "): ";
			test_result |= run_hypervolume_tests(input, output, method_name, eps);
			input.close();
		}
		else {
			output << "Could not open: " << (input_data_testcases_dir + test_name) << std::endl;
		}
		input.close();
	}
	ifs.close();
}

return test_result;
}
