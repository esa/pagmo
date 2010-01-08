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

// 09/09/09 Created by Francesco Biscani.

#include <boost/lexical_cast.hpp>
#include <cstdio>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../../exceptions.h"
#include "../../Functions/rng/rng.h"
#include "../basic/population.h"
#include "base.h"
#include "twodee.h"

namespace pagmo
{
namespace problem {

const char *def_arguments =
    "--controller-type CTRHNN_MULTILAYER "
    "--experiment 28 --sensors D'*' --fitness-function 66816 --evaluation-time-limit 50 "
    "--experiment-arguments sensors=AHE3N,actuators=WG,TauValue=2.5,hiddennodes=10,"
    "disablelogging=yes,arena=9,stepinterval=0.2,maxspeed=8,angle_theshold=1,cameradirectionnoise=0.04,"
    "cameradistancenoise=0.025,cameradirectionbias=0.0,cameradistancebias=0.0,CAM_IR=1,GripperSpeed=1.0,"
    "WheelNoise=0.2 --fitness-function-arguments collisionsallowedpersbot=1 --post-evaluate "
    "--number-of-individuals-to-evaluate 1 --number-of-samples 40 --renderer NULL";

twodee::twodee(int n):base(n),m_random_seed(static_rng_uint32()()),m_arguments(def_arguments) {}

twodee::twodee(int n, const std::string &arguments):base(n),m_random_seed(static_rng_uint32()()),m_arguments(arguments) {}

double twodee::objfun_(const std::vector<double> &v) const
{
	double retval = 0;
	const size_t size = v.size();
	// Write the current chromosome to a temporary file.
	const std::string input_name = std::string(tempnam(0,boost::lexical_cast<std::string>(m_random_seed).c_str()));
	std::ofstream input(input_name.c_str());
	input << size << " ";
	for (size_t i = 0; i < size; ++i) {
		input << v[i];
		if (i + 1 != size) {
			input << " ";
		}
	}
	input.close();
	// Determine the name of the output file.
	const std::string output_name = std::string(tempnam(0,boost::lexical_cast<std::string>(m_random_seed).c_str()));
	// Build and run the command.
	const std::string command = std::string("./twodee ") + std::string("--load-chromosome ") + input_name + std::string(" ") + std::string(m_arguments)
	                            + std::string(" --random-seed ") + boost::lexical_cast<std::string>(m_random_seed)
	                            + std::string(" --twodee_fitness_output ")+ output_name + std::string(" 2>/dev/null");
	int status = std::system(command.c_str());
	std::ifstream res_stream(output_name.c_str());
	if (res_stream.is_open()) {
		std::string line;
		std::getline(res_stream,line);
		retval = boost::lexical_cast<double>(line);
	} else {
		status = 1;
	}
	// Try removing the temporary files.
	std::remove(input_name.c_str());
	std::remove(output_name.c_str());
	if (status) {
		pagmo_throw(std::runtime_error,"error executing twodee");
	}
	return -retval;
}

void twodee::pre_evolution(population & pop) const
{
	m_random_seed = static_rng_uint32()();
	//Re-evaluate the population with respect to the new seed (Internal Sampling Method)
	for (size_t i=0; i<pop.size(); ++i) {
		pop[i] = individual(*this, pop[i].get_decision_vector(), pop[i].get_velocity());
	}
}

}
}
