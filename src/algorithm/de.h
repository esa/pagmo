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

#ifndef PAGMO_ALGORITHM_DE_H
#define PAGMO_ALGORITHM_DE_H

#include <string>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace algorithm {

/// Differential Evolution Algorithm
/**
 *
 * \image html de.jpg "Differential Evolution block diagram."
 * \image latex de.jpg "Differential Evolution block diagram." width=5cm
 *
 * Differential Evolution is an heuristic optimizer developed by Rainer Storn and Kenneth Price.
 *
 * ''A breakthrough happened, when Ken came up with the idea of using vector differences for perturbing
 * the vector population. Since this seminal idea a lively discussion between Ken and Rainer and endless
 * ruminations and computer simulations on both parts yielded many substantial improvements which
 * make DE the versatile and robust tool it is today'' (from the official web pages....)
 *
 * The implementation provided for PaGMO derives from the code provided in the official
 * DE web site and is suitable for box-constrained single-objective continuous optimization.
 *
 * At each call of the evolve method a number of function evaluations equal to m_gen * pop.size()
 * is performed.
 *
 * NOTE: when called on mixed-integer problems DE treats the integer part as fixed and optimizes
 * the continuous part.
 *
 * NOTE2: when called on stochastic optimization problems, DE changes the seed
 * at the end of each generation.
 *
 * NOTE3: the velocity is also updated along DE whenever a new chromosome is accepted.
 *
 * @see http://www.icsi.berkeley.edu/~storn/code.html for the official DE web site
 * @see http://www.springerlink.com/content/x555692233083677/ for the paper that introduces Differential Evolution
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 */

class __PAGMO_VISIBLE de: public base
{
public:
	de(int = 100, double  = 0.8, double = 0.9, int = 2, double = 1e-6, double = 1e-6);
	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
	void set_cr(double cr);
	double get_cr() const;
	void set_f(double cr);
	double get_f() const;
protected:
	std::string human_readable_extra() const;
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<int &>(m_gen);
		ar & const_cast<double &>(m_f);
		ar & const_cast<double &>(m_cr);
		ar & const_cast<double &>(m_ftol);
		ar & const_cast<double &>(m_xtol);
		ar & const_cast<int &>(m_strategy);
	}
	// Number of generations.
	const int m_gen;
	// Weighting factor
	double m_f;
	// Crossover probability
	double m_cr;
	// Startegy
	const int m_strategy;
	const double m_ftol;
	const double m_xtol;
	
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::de)

#endif // DE_H
