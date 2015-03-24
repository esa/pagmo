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

#ifndef PAGMO_ALGORITHM_NSGA2_H
#define PAGMO_ALGORITHM_NSGA2_H

#include "../config.h"
#include "../serialization.h"
#include "base.h"


namespace pagmo { namespace algorithm {

/// Nondominated Sorting genetic algorithm II (NSGA-II)
/**
 * NSGA-II is a nondominated-sorting based multiobjective evolutionary algorithm.
 * It genererates offspring with crossover and mutation and select the next
 * generation according to nondominated-sorting and crowding distance comparison.
 *
 * The algorithm can be applied to continuous box-bounded optimization. The version for mixed integer
 * and constrained optimization is also planned.
 * 
 * @see Deb, K. and Pratap, A. and Agarwal, S. and Meyarivan, T., "A fast and elitist multiobjective genetic algorithm: NSGA-II"
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 */

class __PAGMO_VISIBLE nsga2: public base
{
public:
	nsga2(int gen=100, double cr = 0.95, double eta_c = 10, double m = 0.01, double eta_m = 50);
	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
	
protected:
	std::string human_readable_extra() const;
	
private:
	pagmo::population::size_type tournament_selection(pagmo::population::size_type, pagmo::population::size_type, const pagmo::population&) const;
	void crossover(decision_vector&, decision_vector&, pagmo::population::size_type, pagmo::population::size_type,const pagmo::population&) const;
	void mutate(decision_vector&, const pagmo::population&) const;
	
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<int &>(m_gen);
		ar & const_cast<double &>(m_cr);
		ar & const_cast<double &>(m_eta_c);
		ar & const_cast<double &>(m_m);
		ar & const_cast<double &>(m_eta_m);
	}
	//Number of generations
	const int m_gen;
	//Crossover rate
	const double m_cr;
	// Ditribution index for crossover
	const double m_eta_c;
	// Mutation rate
	const double m_m;
	// Ditribution index for mutation
	const double m_eta_m;


};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::nsga2)

#endif // PAGMO_ALGORITHM_NSGA2_H
