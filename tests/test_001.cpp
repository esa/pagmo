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
#include <boost/iterator/transform_iterator.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/utility/result_of.hpp>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

#include "../src/algorithm/de.h"
#include "../src/island.h"
#include "test_functions.h"

using namespace pagmo;

struct std_calculator {
	std_calculator(const double &mean):m_mean(mean) {}
	typedef double result_type;
	double operator()(const double &x) const
	{
		return (m_mean - x) * (m_mean - x);
	}
	const double m_mean;
};

typedef boost::transform_iterator<std_calculator,std::vector<double>::iterator> std_iterator;

int main()
{
	algorithm::de de = algorithm::de(500,.8,.8,2);
	const std::vector<problem::base_ptr> probs(get_test_problems());
	std::cout << "Testing algorithm: " << de << '\n';
	for (std::vector<problem::base_ptr>::const_iterator it = probs.begin(); it != probs.end(); ++it) {
		std::cout << "\tTesting problem: " << (**it) << '\n';
		std::vector<double> champs;
		for (int i = 0; i < 100; ++i) {
			island isl(**it,de,20,1);
			isl.evolve(1);
			isl.join();
			champs.push_back(isl.get_population().champion().f[0]);
		}
		std::cout << "\t\tBest:\t" << boost::lexical_cast<std::string>(*std::min_element(champs.begin(),champs.end())) << '\n';
		const double mean = std::accumulate(champs.begin(),champs.end(),double(0.)) / champs.size();
		std::cout << "\t\tMean:\t" << boost::lexical_cast<std::string>(mean) << '\n';
		std::cout << "\t\tStd:\t" <<  boost::lexical_cast<std::string>(std::sqrt(
			std::accumulate(
				std_iterator(champs.begin(),std_calculator(mean)),
				std_iterator(champs.end(),std_calculator(mean)),
				double(0.)
			) / champs.size())
		) << '\n';
	}
	return 0;
}
