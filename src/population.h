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

#ifndef PAGMO_POPULATION_H
#define PAGMO_POPULATION_H

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "config.h"
#include "problem/base.h"
#include "rng.h"
#include "types.h"

namespace pagmo
{

// Forward declaration of island class, needed for friendship.
class island;

/// Population class.
/**
 * This class contains an instance of an optimisation problem and a group of candidate solutions represented by the class individual_type. On creation,
 * the population is associated to a problem and initialised with random decision vectors. An instance of champion_type automatically keeps track of the best solution
 * ever appeared in the population. Methods are offered to get and manipulate the single individuals.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE population
{
		friend class island;
	public:
		/// Individuals stored in the population.
		/**
		 * Individuals store the current decision and velocity vectors, the current constraint vector and the current fitness vector. They also
		 * keep memory of the best decision, constraint and fitness vectors "experienced" so far by the individual.
		 */
		struct individual_type {
			/// Current decision vector.
			decision_vector		cur_x;
			/// Current velocity vector.
			decision_vector		cur_v;
			/// Current constraint vector.
			constraint_vector	cur_c;
			/// Current fitness vector.
			fitness_vector		cur_f;
			/// Best decision vector so far.
			decision_vector		best_x;
			/// Best constraint vector so far.
			constraint_vector	best_c;
			/// Best fitness vector so far.
			fitness_vector		best_f;
		};
		/// Population champion.
		/**
		 * A champion is the best individual that ever lived in the population. It is defined by a decision vector, a constraint vector and a fitness vector.
		 */
		struct champion_type {
			/// Decision vector.
			decision_vector		x;
			/// Constraint vector.
			constraint_vector	c;
			/// Fitness vector.
			fitness_vector		f;
		};
		/// Population size type.
		typedef std::vector<individual_type>::size_type size_type;
		/// Const iterator.
		typedef std::vector<individual_type>::const_iterator const_iterator;
		population(const problem::base &, int n = 0);
		population(const population &);
		population &operator=(const population &);
		const individual_type &get_individual(const size_type &) const;
		const problem::base &problem() const;
		const champion_type &champion() const;
		std::string human_readable_terse() const;
		std::string human_readable() const;
		size_type get_best_idx() const;
		size_type get_worst_idx() const;
		void set_x(const size_type &, const decision_vector &);
		void set_v(const size_type &, const decision_vector &);
		size_type size() const;
		const_iterator begin() const;
		const_iterator end() const;
		/// Compare individuals according to their current fitness and constraints.
		/**
		 * This functor is used to sort the individuals according to their current fitness-constraint vector pair.
		 *
		 * @see problem::base::compare_fc.
		 */
		struct cur_fc_comp {
			/// Constructor from problem.
			/**
			 * p.compare_fc() will be used to rank the individuals.
			 *
			 * @param[in] p problem::base which will provide the ranking method.
			 */
			cur_fc_comp(const problem::base &p):m_p(p) {}
			/// Compare individuals.
			/**
			 * @param[in] i1 first individual.
			 * @param[in] i2 second individual.
			 *
			 * @return true if i1 is strictly better than i2, false otherwise.
			 */
			bool operator()(const population::individual_type &i1, const population::individual_type &i2) const
			{
				return m_p.compare_fc(i1.cur_f,i1.cur_c,i2.cur_f,i2.cur_c);
			}
			/// Const reference to the problem that provides the ranking method.
			const problem::base &m_p;
		};
	private:
		population();
	private:
		typedef std::vector<individual_type> container_type;
		// Data members.
		// Problem.
		problem::base_ptr	m_prob;
		// Container of individuals.
		container_type		m_container;
		// Population champion.
		champion_type		m_champion;
		// Double precision random number generator.
		mutable	rng_double	m_drng;
		// uint32 random number generator.
		mutable	rng_uint32	m_urng;
};

__PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &, const population &);

}

#endif
