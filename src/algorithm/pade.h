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

#ifndef PAGMO_ALGORITHM_PADE_H
#define PAGMO_ALGORITHM_PADE_H

#include "../config.h"
#include "../serialization.h"
#include "base.h"
#include "jde.h"
#include "../problem/decompose.h"



namespace pagmo { namespace algorithm {

/// Parallel Decomposition (PaDe)
/**
 *
 * This class implement a multi-objective optimization algorithm based on parallel decomposition.
 * For each element of the population a different single objective problem is generated using
 * a decomposition method. Those single-objective problems are thus solved in parallel.
 * At the end of the evolution the population is set as the best individual for each single-objective problem.
 *
 * PaDe assumes all the objectives need to be minimized.
 *
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 * @author Dario Izzo (dario.izzo@gmail.com)
 **/

class __PAGMO_VISIBLE pade: public base
{
public:
	/// Mechanism used to generate the weight vectors
	enum weight_generation_type {
	    RANDOM=0, ///< Weights are generated uniformly at random on the simplex 
	    GRID=1,///< Weights are generated on a uniform grid layed down on the simplex
	    LOW_DISCREPANCY=2 ///< Weights are generated on the simplex with low-discrepancy
	};
	pade(
		 int gen=10, 
		 unsigned int max_parallelism = 1, 
		 pagmo::problem::decompose::method_type =  pagmo::problem::decompose::BI,
		 const pagmo::algorithm::base & = pagmo::algorithm::jde(100),
		 population::size_type = 8,
		 weight_generation_type = LOW_DISCREPANCY, 
		 const fitness_vector & = std::vector<double>()
		);
	pade(const pade &);

	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
	std::vector<fitness_vector> generate_weights(const unsigned int, const unsigned int) const;

protected:
	std::string human_readable_extra() const;

private:
	void reksum(std::vector<std::vector<double> > &, const std::vector<unsigned int>&, unsigned int, unsigned int, std::vector<double> = std::vector<double>() ) const;
	void compute_neighbours(std::vector<std::vector<int> > &, const std::vector<std::vector <double> > &);
	double distance(pagmo::fitness_vector , pagmo::fitness_vector);
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<int &>(m_gen);
		ar & const_cast<unsigned int &>(m_threads);
		ar & const_cast<pagmo::problem::decompose::method_type &>(m_method);
		ar & const_cast<base_ptr &>(m_solver);
		ar & const_cast<population::size_type &>(m_T);
		ar & const_cast<weight_generation_type &>(m_weight_generation);
		ar & m_z;
	}
	//Number of generations
	const int m_gen;
	const unsigned int m_threads;
	const pagmo::problem::decompose::method_type m_method;
	const base_ptr m_solver;
	const population::size_type m_T;
	const weight_generation_type m_weight_generation;
	fitness_vector m_z;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::pade)

#endif // PAGMO_ALGORITHM_PADE_H
