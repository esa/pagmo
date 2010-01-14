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

// 04/06/08 Created by Dario Izzo.

#ifndef PAGMO_PROBLEM_BASE_H
#define PAGMO_PROBLEM_BASE_H

#include <iostream>
#include <string>
#include <typeinfo>
#include <vector>

#include "../../config.h"
#include "../../atomic_counters/atomic_counters.h"

namespace pagmo
{

// Forward declaration of population class. We cannot include the header directly because
// we would incur into circular dependency problem (population.h includes base.h which
// includes population.h ...)
class population;

/// Problem namespace.
/**
 * This namespace contains all the problems implemented in PaGMO.
 **/
namespace problem {

/// Base problem class.
/**
 * This class implements the highest hierarchical level of a global optimisation problem here defined as unconstrained
 * optimisation problems where the search domain is hyperrectangular.
 */
class __PAGMO_VISIBLE base
{
		friend std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const base &);
		friend size_t __PAGMO_VISIBLE_FUNC objfun_calls();
		friend void __PAGMO_VISIBLE_FUNC reset_objfun_calls();
	public:
		// Bounds getters and setters via reference
		const std::vector<double> &get_lb() const;
		const std::vector<double> &get_ub() const;
		void set_lb(const std::vector<double> &);
		void set_lb(size_t, const double &);
		void set_lb(size_t, const double *);
		void set_ub(const std::vector<double> &);
		void set_ub(size_t, const double &);
		void set_ub(size_t, const double *);		
		// Dimension getter
		size_t getDimension() const;
		double objfun(const std::vector<double> &) const;
		virtual base *clone() const = 0;
		std::string id_name() const;
		virtual void pre_evolution(population &) const {}
		virtual void post_evolution(population &) const {}
		virtual bool operator==(const base &) const;
		bool operator!=(const base &) const;

		/// Get the name identyfing the object (<b>not</b> the class).
		/** Exposed to Python. The string should identify the object, so that instantiations of the same class with different parameters are distinguishable. */
		virtual std::string id_object() const = 0;

	protected:
		/// Print function.
		/**
		 * Called by operator<<, can be re-implemented in subclasses. Default implementation prints the problem's C++ name,
		 * dimension and bounds.
		 */
		virtual std::ostream &print(std::ostream &) const;
		// The objective function - must be implemented in subclasses
		virtual double objfun_(const std::vector<double> &) const = 0;
		// Constructor from size; construct problem of size n, with lower bounds to zero and upper bounds to one
		base(int);
		// Constructor with array bounds initialisers
		base(const size_t &, const double *, const double *);
		// Constructor with vectors initialisers
		base(const std::vector<double> &, const std::vector<double> &);
		// These need to be protected and cannot be private and const as some problems need to deifne the LB and UB at run time (i.e. LJ)
		std::vector<double> LB;
		std::vector<double> UB;
	private:
		void check_boundaries() const;
		static atomic_counter_size_t m_objfun_counter;
};

std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const base &);

size_t __PAGMO_VISIBLE_FUNC objfun_calls();
void __PAGMO_VISIBLE_FUNC reset_objfun_calls();

}
}

#endif
