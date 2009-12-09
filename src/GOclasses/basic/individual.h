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

// 16/05/08 Created by Dario Izzo.

#ifndef PAGMO_INDIVIDUAL_H
#define PAGMO_INDIVIDUAL_H

#include <iostream>
#include <vector>

#include "../../../config.h"
#include "../../Functions/rng/rng.h"
#include "../../exceptions.h"
#include "../problems/GOproblem.h"

/// Individual class.
/**
 * Individuals represent solutions for a problem.
 * \todo The desturbing thing is passing the problem as an argument of the constructor. It suggests, that
 * the Individual keeps an internal reference to it, which is dangerous (GOProblem is not assumed to be thread safe).
 * It the end the reference is not kept, but in order to discover this, one must re-read the code of the
 * Individual.
 */
class __PAGMO_VISIBLE Individual {
public:
	/// Constructor.
	/**
	 * Constructs an individual for the given problem with x randomly placed within the problem bounds
	 * and a random velocity of maximum magnitude (UB-LB).
	 * \param[in] problem concerned problem.
	 */
	Individual(const GOProblem& problem);
	
	/// Constructor.
	/**
	 * Constructs an Individual with given x and v.
	 * \param[in] problem concerned problem.
	 * \param[in] x_ individual position.
	 * \param[in] v_ individual velocity.
	 */
	Individual(const GOProblem& problem, const std::vector<double>& x_, const std::vector<double>& v_);
	
	/// Constructor.
	/**
	 * Constructs an Individual with given x and an empty velocity.
	 * \param[in] problem concerned problem.
	 * \param[in] x_ individual position.
	 */
	Individual(const GOProblem& problem, const std::vector<double> & x_);
	
	/// Constructor.
	/**
	 * Constructs an Individual with given position, velocity and fitness.
	 * \param[in] x_ individual position.
	 * \param[in] v_ individual velocity.
	 * \param[in] fitness_ individual fitness.
	 */
	Individual(const std::vector<double> &x_, const std::vector<double> &v_, const double &fitness_)
			:x(x_),
			v(v_),
			fitness(fitness_)
	{
		if (x.size() != v.size()) {
			pagmo_throw(value_error,"while constructing individual, size mismatch between decision vector and velocity vector");
		}
	}
	
	/// Constructor.
	/**
	 * Constructs an Individual with given position and fitness.
	 * \param[in] x_ individual position.
	 * \param[in] fitness_ individual fitness.
	 */
	Individual(const std::vector<double> &x_, const double &fitness_)
			:x(x_),
			v(x_.size()),
			fitness(fitness_)
	{
	}
	
	/// Copy constructor.
	/**
	 * Creates a deep copy of an individual.
	 * \param[in] individual The individual to be duplicated.
	 */
	Individual(const Individual& individual)
		:x(individual.x),
		v(individual.v),
		fitness(individual.fitness)
	{
	}
	
	/// Assignment operator.
	Individual& operator=(const Individual &);
	
	///Returns the Individual fitness.
	double getFitness() const { return fitness; }
	
	///Returns the Individual chromosome (position).
	const std::vector<double> &getDecisionVector() const { return x; }
	
	///Returns the Individual velocity.
	const std::vector<double> &getVelocity() const { return v; }
	
	/// Check if the individual is compatible with a given problem.
	/**
	 * Compatibility means that the individual's size is the same as the problem's and that the values
	 * of the decision vector are within the boundaries defined in the problem. If the individual is incompatible,
	 * an exception will be thrown.
	 * \param p Problem of interest.
	 */
	void check(const GOProblem &p) const;
	
	/// Individual comparator for STL sorting methods.
	/**
	 * Useful utility function. Minimisation is assumed (lower fitness goes first).
	 */
	static int compare_by_fitness(const Individual& ind1, const Individual& ind2);
	
private:
	/// Stream output operator
	friend std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const Individual &);
	
	std::vector<double>			x; ///< Individual chromosome (position).
	std::vector<double>			v; ///< Individual velocity.
	double						fitness; ///< Individual fitness.
};

/// Stream output operator
std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const Individual &);

#endif
