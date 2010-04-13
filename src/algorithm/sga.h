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

#ifndef PAGMO_ALGORITHM_SGA_H
#define PAGMO_ALGORITHM_SGA_H

#include "../config.h"
#include "base.h"
#include "../problem/base.h"


namespace pagmo { namespace algorithm {

/// The Simple Genetic Algorithm (SGA)
/**
 * Genetic algorithms are very popular algorithms used widely by people of very different backgrounds.
 * As a consequence there are a large number of different implementations and toolboxes that are available
 * and can be used to construct a genetic algorithm. We decided not to choose one of these and, instead, to
 * provide only a basic implementation of the algorithm implementing a floating point encoding (not binary)
 * and some common mutation and crossover strategies, hence the name Simple Genetic Algorithm.
 *
 * Mutation is gaussian or random, crossover exponential or binomial and selection is tournament or
 * best20 (i.e. 20% best of the population is selcted and reproduced 5 times).
 *
 * The algorithm works on single objective, box constrained problems. The mutation operator acts
 * differently on continuous and discrete variables.
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 *
 */

class __PAGMO_VISIBLE sga: public base
{
public:
	/// Selection info
	struct selection{ enum type {BEST20 = 0,ROULETTE = 1}; };
	/// Mutation operator info
	struct mutation {
		enum type {GAUSSIAN = 0, RANDOM = 1};
		mutation(mutation::type t, double width) : m_type(t),m_width(width) {}
		type m_type;
		double m_width;
	};

	/// Crossover operator info
	struct crossover { enum type {BINOMIAL = 0, EXPONENTIAL = 1}; };
	sga(int gen, const double &cr, const double &m, int elitism = 1,
	    mutation::type mut  = mutation::GAUSSIAN, double width = 0.05,
	    selection::type sel = selection::ROULETTE,
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
	const double& m_cr;
	//Mutation rate
	const double& m_m;
	//Elitism (number of generations after which to reinsert the best)
	const int m_elitism;
	//Mutation
	const mutation m_mut;
	//Selection_type
	const selection::type m_sel;
	//Crossover_type
	const crossover::type m_cro;


};

}} //namespaces

#endif // PAGMO_ALGORITHM_SGA_H
