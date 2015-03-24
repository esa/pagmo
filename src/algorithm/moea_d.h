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

#ifndef PAGMO_ALGORITHM_MOEAD_H
#define PAGMO_ALGORITHM_MOEAD_H

#include "../config.h"
#include "../serialization.h"
#include "base.h"
#include "../problem/decompose.h"



namespace pagmo { namespace algorithm {

/// MOEA/D - DE
/**
 *
 * This class implements the Multi Objective Evolutionary Algorithm based on Decomposition and Differential Evolution
 * crossover. The reference papers can be found below. By activating or deactivating the bool preserve_diversity
 * one can select to use the ideas introduced in the second paper or not. In all cases Tchebycheff decomposition and 
 * a differential evolution operator are used
 *
 * @see Zhang, Qingfu, and Hui Li. "MOEA/D: A multiobjective evolutionary algorithm based on decomposition." Evolutionary Computation, IEEE Transactions on 11.6 (2007): 712-731.
 * @see Li, Hui, and Qingfu Zhang. "Multiobjective optimization problems with complicated Pareto sets, MOEA/D and NSGA-II." Evolutionary Computation, IEEE Transactions on 13.2 (2009): 284-302.
 *
 **/

class __PAGMO_VISIBLE moead: public base
{
public:
	/// Mechanism used to generate the weight vectors
	enum weight_generation_type {
	 RANDOM=0, ///< Weights are generated uniformly at random on the simplex 
	 GRID=1,///< Weights are generated on a uniform grid layed down on the simplex
	 LOW_DISCREPANCY=2 ///< Weights are generated on the simplex with low-discrepancy
	};
	
	moead(
		 int gen=100, 
		 weight_generation_type = GRID, 
		 population::size_type = 20,
		 double realb = 0.9,
		 unsigned int limit = 2,
		 double CR = 1.0,
		 double F=0.5,
		 double eta_m = 20,
		 bool preserve_diversity = true
		);

	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
	std::vector<fitness_vector> generate_weights(const unsigned int, const unsigned int) const;

protected:
	std::string human_readable_extra() const;

private:
	void reksum(std::vector<std::vector<double> > &, const std::vector<unsigned int>&, unsigned int, unsigned int, std::vector<double> = std::vector<double>() ) const;
	void compute_neighbours(std::vector<std::vector<int> > &, const std::vector<std::vector <double> > &);
	void mating_selection(std::vector<population::size_type> &, int, int,const std::vector<std::vector<population::size_type> >&) const;
	void mutation(decision_vector&, const population&, double rate) const;
	
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<int &>(m_gen);
		ar & const_cast<population::size_type &>(m_T);
		ar & const_cast<weight_generation_type &>(m_weight_generation);
		ar & const_cast<double &>(m_realb);
		ar & const_cast<unsigned int &>(m_limit);
		ar & const_cast<double &>(m_cr);
		ar & const_cast<double &>(m_f);
		ar & const_cast<double &>(m_eta_m);
		ar & const_cast<double &>(m_preserve_diversity);
	}
	//Number of generations
	const int m_gen;
	const population::size_type m_T;
	const weight_generation_type m_weight_generation;
	//probability of selecting mating parents from neighborhood
	const double m_realb;
	const unsigned int m_limit;
	const double m_cr;
	const double m_f;
	const double m_eta_m;
	const double m_preserve_diversity;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::moead)

#endif // PAGMO_ALGORITHM_MOEAD_H
