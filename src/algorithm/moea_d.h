/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
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

#ifndef PAGMO_ALGORITHM_MOEAD_H
#define PAGMO_ALGORITHM_PADE_H

#include "../config.h"
#include "../serialization.h"
#include "base.h"
#include "../problem/decompose.h"



namespace pagmo { namespace algorithm {

/// MOEA/D
/**
 *
 * This class implements the 
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
		 pagmo::problem::decompose::method_type =  pagmo::problem::decompose::TCHEBYCHEFF,
		 population::size_type = 20,
		 weight_generation_type = GRID, 
		 double realb = 0.9,
		 unsigned int limit = 2
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
		ar & const_cast<pagmo::problem::decompose::method_type &>(m_method);
		ar & const_cast<population::size_type &>(m_T);
		ar & const_cast<weight_generation_type &>(m_weight_generation);
	}
	//Number of generations
	const int m_gen;
	const pagmo::problem::decompose::method_type m_method;
	const population::size_type m_T;
	const weight_generation_type m_weight_generation;
	//probability of selecting mating parents from neighborhood
	const double m_realb;
	const unsigned int m_limit;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::moead)

#endif // PAGMO_ALGORITHM_MOEAD_H
