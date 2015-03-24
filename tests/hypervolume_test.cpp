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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <locale>
#include <boost/date_time/posix_time/posix_time.hpp>
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

/// Used for generating the report files and outputing to the std::cout.
class dual_stream
{
public:
	dual_stream(std::ostream &os1, std::ostream &os2) : m_os1(os1), m_os2(os2) {
		m_os1.precision(15);
		m_os2.precision(15);
	}
	template<class T>
	dual_stream& operator<<(const T& item) {
		m_os1 << item;
		m_os2 << item;
		return *this;
	}
private:
	std::ostream &m_os1;
	std::ostream &m_os2;

};

class hypervolume_test
{
public:
	hypervolume_test(std::istream &input, dual_stream &output, std::string test_type, std::string method_name, double eps)
		: m_input(input), m_output(output), m_test_type(test_type), m_eps(eps) {

		// create correct algorithm object
		if (method_name == "hv2d") {
			m_method = util::hv_algorithm::hv2d().clone();
		} else if (method_name == "hv3d") {
			m_method = util::hv_algorithm::hv3d().clone();
		} else if (method_name == "hv4d") {
			m_method = util::hv_algorithm::hv4d().clone();
		} else if (method_name == "wfg") {
			m_method = util::hv_algorithm::wfg().clone();
		} else if (method_name == "fpl") {
			m_method = util::hv_algorithm::fpl().clone();
		} else if (method_name == "hoy") {
			m_method = util::hv_algorithm::hoy().clone();
		} else if (method_name == "bf_approx") {
			m_method = util::hv_algorithm::bf_approx().clone();
		} else if (method_name == "bf_fpras") {
			m_method = util::hv_algorithm::bf_fpras(0.1, 0.1).clone();
		} else {
			output << "Unknown method (" << method_name << ") .. exiting\n";
			exit(1);
		}
	}

	int run_test() {
		m_input >> m_num_tests;

		int OK_counter=0;
		for(int t=0 ; t < m_num_tests ; ++t) {
			load_common();
			util::hypervolume hv_obj = util::hypervolume(m_points);

			// run correct test
			if (m_test_type == "compute") {
				load_compute();
				double hypvol = hv_obj.compute(m_ref_point, m_method);
				if (fabs(hypvol-m_hv_ans) < m_eps) {
					++OK_counter;
				} else {
					m_output << "\n Error in test " << t << ". Got: " << hypvol << ", Expected: " << m_hv_ans << " (abs error: " << fabs(hypvol-m_hv_ans) << ")";
				}

			} else if (m_test_type == "exclusive") {
				load_exclusive();
				double hypvol = hv_obj.exclusive(m_p_idx, m_ref_point, m_method);
				if (fabs(hypvol-m_hv_ans) < m_eps) {
					++OK_counter;
				} else {
					m_output << "\n Error in test " << t << ". Got: " << hypvol << ", Expected: " << m_hv_ans << " (abs error: " << fabs(hypvol-m_hv_ans) << ")";
				}

			} else if (m_test_type == "least_contributor") {
				load_least_contributor();
				unsigned int point_idx = hv_obj.least_contributor(m_ref_point, m_method);
				if (point_idx == m_idx_ans) {
					++OK_counter;
				} else {
					m_output << "\n Error in test " << t << ". Got: " << point_idx << ", Expected: " << m_idx_ans ;
				}

			} else if (m_test_type == "greatest_contributor") {
				load_least_contributor(); // loads the same data as least contributor
				unsigned int point_idx = hv_obj.greatest_contributor(m_ref_point, m_method);
				if (point_idx == m_idx_ans) {
					++OK_counter;
				} else {
					m_output << "\n Error in test " << t << ". Got: " << point_idx << ", Expected: " << m_idx_ans;
				}

			} else {
				m_output << "Unknown test type (" << m_test_type << ") .. exiting";
				exit(1);
			}
		}
		m_output << "\n" << " " << OK_counter << "/" << m_num_tests << " passed";
		return (OK_counter < m_num_tests ? 1 : 0);
	}
private:
	void load_common()
	{
		m_input >> m_f_dim >> m_num_points;
		m_points = std::vector<fitness_vector>(m_num_points, fitness_vector(m_f_dim, 0.0));
		m_ref_point = fitness_vector(m_f_dim, 0.0);
		for(int d = 0 ; d < m_f_dim ; ++d) {
			m_input >> m_ref_point[d];
		}
		for(int i = 0 ; i < m_num_points ; ++i) {
			for(int d = 0 ; d < m_f_dim ; ++d) {
				m_input >> m_points[i][d];
			}
		}
	}

	void load_compute()
	{
		m_input >> m_hv_ans;
	}

	void load_exclusive()
	{
		m_input >> m_p_idx;
		m_input >> m_hv_ans;
	}

	void load_least_contributor()
	{
		m_input >> m_idx_ans;
	}

	int m_num_tests, m_f_dim, m_num_points, m_p_idx;
	unsigned int m_idx_ans;
	double m_hv_ans;
	fitness_vector m_ref_point;
	std::vector<fitness_vector> m_points;

	util::hv_algorithm::base_ptr m_method;
	std::istream &m_input;
	dual_stream &m_output;
	std::string m_test_type;
	double m_eps;
};

int main(int argc, char *argv[])
{
	std::string line;

	// root directory of the hypervolume data
	std::string input_data_dir("hypervolume_test_data/");
	// root directory of the testcases
	std::string input_data_testcases_dir(input_data_dir + "testcases/");
	// testcases list filename
	std::ifstream ifs;
	if (argc > 1) {
		std::string testcases_file = argv[1];
		ifs.open(testcases_file.c_str());
	} else {
		ifs.open((input_data_dir + "testcases_list.txt").c_str());
	}
	std::ofstream report_stream("hypervolume_test_report.txt");
	dual_stream output(report_stream, std::cout);

	// date and time for the reference
	output << boost::posix_time::second_clock::local_time() << "\n";
	int test_result = 0;
	if (ifs.is_open()) {
		while (ifs.good()) {
			std::string line;
			getline(ifs,line);
			if (line == "" || line[0] == '#')
				continue;
			std::stringstream ss (line);
			std::string test_type;
			std::string method_name;
			std::string test_name;
			double eps;
			ss >> test_type;
			ss >> method_name;
			ss >> test_name;
			ss >> eps;
			std::ifstream input((input_data_testcases_dir + test_name).c_str());
			if (input.is_open()){
				output << test_type << " / " << method_name << " / " << test_name << " / eps:" << eps << "\n";

				hypervolume_test hvt(input, output, test_type, method_name, eps);

				boost::posix_time::ptime time_start(boost::posix_time::microsec_clock::local_time());
				test_result |= hvt.run_test();
				boost::posix_time::ptime time_end(boost::posix_time::microsec_clock::local_time());
				boost::posix_time::time_duration time_diff(time_end - time_start);

				output << " (Time " << time_diff.total_milliseconds() / 1000.0 << " s)\n\n";

				input.close();
			}
			else {
				output << "Could not open: " << (input_data_testcases_dir + test_name) << "\n";
			}
			input.close();
		}
		ifs.close();
	}

	return test_result;
}
