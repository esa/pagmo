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

#ifndef PAGMO_ALGORITHM_SGA_GRAY_H
#define PAGMO_ALGORITHM_SGA_GRAY_H

#include "../config.h"
#include "../problem/base.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace algorithm {

/// The Simple Genetic Algorithm with gray encoding (sga_gray)
/**
 * Simple Genetic algorithms are very popular algorithms used widely by people of very
 * different backgrounds.
 * Contrary to PAGMO::SGA, this implementation uses a binary encoding and some common
 * mutation and crossover strategies.
 *
 * Mutation is random, crossover uniform and selection is roulette or
 * best20 (i.e. 20% best of the population is selected and reproduced 5 times).
 *
 * The algorithm works on single objective, box constrained problems.
 *
 * @author Jeremie Labroquere (jeremie.labroquere@gmail.com)
 */

class __PAGMO_VISIBLE sga_gray: public base
{
public:
	/// Selection info
	struct selection {
		/// Selection type, best 20% or roulette
		enum type {BEST20 = 0,ROULETTE = 1};
	};
	/// Mutation operator info
	struct mutation {
		/// Mutation type, uniform
		enum type {UNIFORM = 0};
	};

	/// Crossover operator info
	struct crossover {
		/// Crossover type, single point
		enum type {SINGLE_POINT = 0};
	};
	sga_gray(int gen  = 1, const double &cr = .95, const double &m = .02, int elitism = 1,
			 mutation::type mut  = mutation::UNIFORM,
			 selection::type sel = selection::ROULETTE,
			 crossover::type cro = crossover::SINGLE_POINT);
	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
protected:
	std::string human_readable_extra() const;
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<int &>(m_gen);
		ar & const_cast<double &>(m_cr);
		ar & const_cast<double &>(m_m);
		ar & const_cast<int &>(m_elitism);
		ar & const_cast<mutation::type &>(m_mut);
		ar & const_cast<selection::type &>(m_sel);
		ar & const_cast<crossover::type &>(m_cro);
		ar & m_max_encoding_integer;
		ar & m_bit_encoding;
	}
	//Number of generations
	const int m_gen;
	//Crossover rate
	const double m_cr;
	//Mutation rate
	const double m_m;

	//Elitism (number of generations after which to reinsert the best)
	const int m_elitism;
	//Mutation
	const mutation::type m_mut;
	//Selection_type
	const selection::type m_sel;
	//Crossover_type
	const crossover::type m_cro;

private:
	// genetic algoritms operators
	std::vector<int> selection(const std::vector<fitness_vector> &, const problem::base &) const;
	void crossover(std::vector< std::vector<int> > &pop_x) const;
	void mutate(std::vector< std::vector<int> > &pop_x) const;

private:
	// binary conversion
	std::vector<int> double_to_binary(const double &number, const double &lb, const double &ub) const;
	double binary_to_double(const std::vector<int> &binary, const double &lb, const double &ub) const;
	std::vector<int> gray_to_binary(const std::vector<int> &gray) const;
	std::vector<int> binary_to_gray(const std::vector<int> &binary) const;

	// encoding/decoding
	std::vector<int> encode(const decision_vector &x_vector, const decision_vector &lb, const decision_vector &ub) const;
	decision_vector decode(const std::vector<int> &chrom, const decision_vector &lb, const decision_vector &ub) const;

	// encoding size
	int m_max_encoding_integer;
	int m_bit_encoding;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::sga_gray)

#endif // PAGMO_ALGORITHM_SGA_GRAY_H
