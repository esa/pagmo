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

#ifndef PAGMO_UTIL_RACE_ALGO_H
#define PAGMO_UTIL_RACE_ALGO_H

#include <iostream>
#include <string>
#include <vector>

#include "../config.h"
#include "../serialization.h"
#include "../problem/base.h"
#include "../problem/ackley.h"
#include "../algorithm/base.h"

namespace pagmo { namespace util { namespace racing {

/// Racing mechanism for algorithms
/**
 * This class allows the racing of a set of algorithms on a problem or a set of
 * problems. It supports the racing over single objective box-constrained and
 * equality / inequality constrained problems.
 */
class __PAGMO_VISIBLE race_algo
{
	public:
		race_algo(const std::vector<algorithm::base_ptr> &algos = std::vector<algorithm::base_ptr>(), const problem::base &prob = problem::ackley(), unsigned int pop_size = 100, unsigned int seed = 0);
		race_algo(const std::vector<algorithm::base_ptr> &algos, const std::vector<problem::base_ptr> &prob, unsigned int pop_size = 100, unsigned int seed = 0);

		// Main method containing all the juice
		std::pair<std::vector<unsigned int>, unsigned int> run(
			const unsigned int n_final,
			const unsigned int min_trials,
			const unsigned int max_count,
			double delta,
			const std::vector<unsigned int> &,
			const bool race_best,
			const bool screen_output
		);

	private:

		std::vector<algorithm::base_ptr> m_algos;
		std::vector<problem::base_ptr> m_probs;
		unsigned int m_pop_size;
		unsigned int m_seed;
};

}}}

#endif
