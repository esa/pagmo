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

#ifndef PAGMO_ALGORITHM_NSGA2_H
#define PAGMO_ALGORITHM_NSGA2_H

#include "../config.h"
#include "base.h"
#include "../problem/base.h"


namespace pagmo { namespace algorithm {

/// Nondominated Sorting genetic algorithm II (NSGA-II)
/**
 * NSGA-II is a nondominated-sorting based multiobjective evolutionary algorithm.
 * It genererates offspring with crossover and mutation and select the next
 * generation according to nondominated-sorting and crowding distance comparison.
 *
 * The algorithm works on continuous box-constrained multi objective problems.
 *
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 */

class __PAGMO_VISIBLE nsga2: public base
{
public:
	/// Mutation operator info
	struct mutation {
		/// Mutation type, gaussian or random
		enum type {GAUSSIAN = 0, RANDOM = 1};
		/// Constructor
		/**
		 * \param[in] t the mutation type
		 * \param[in] width the width of the gaussian bell in case of a gaussian mutation. The
		 *		parameter is otherwise ignored. width is a percentage with respect to the
		 *		ub[i]-lb[i] width.
		 */
		mutation(mutation::type t, double width) : m_type(t),m_width(width) {}
		/// Mutation type
		type m_type;
		/// Mutation width
		double m_width;
	};

	/// Crossover operator info
	struct crossover {
		/// Crossover type, binomial or exponential
		enum type {BINOMIAL = 0, EXPONENTIAL = 1};
	};

	nsga2(int gen, const double &cr, const double &m,
	    mutation::type mut  = mutation::GAUSSIAN, double width = 0.05,
	    crossover::type cro = crossover::EXPONENTIAL);
	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
protected:
	std::string human_readable_extra() const;
private:
	//Number of generations
	const int m_gen;
	//Crossover rate
	const double m_cr;
	//Mutation rate
		const double m_m;
	//Mutation
	const mutation m_mut;
	//Crossover_type
	const crossover::type m_cro;


};

}} //namespaces

#endif // PAGMO_ALGORITHM_NSGA2_H
