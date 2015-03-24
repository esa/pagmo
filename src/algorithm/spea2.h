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

#ifndef PAGMO_ALGORITHM_SPEA2_H
#define PAGMO_ALGORITHM_SPEA2_H

#include "../config.h"
#include "../serialization.h"
#include "../util/neighbourhood.h"
#include "base.h"

namespace pagmo { namespace algorithm {

class  distance_sorter {
public:
	distance_sorter(const std::vector<std::vector<pagmo::population::size_type>  > &neighbours,
					const std::vector<fitness_vector> &fit):
		m_neighbours(neighbours),
		m_fit(fit) {}
	bool operator()(unsigned int a, unsigned int b) {
		if(a>=m_fit.size() || a>=m_neighbours.size()){
			pagmo_throw(value_error,"SPEA2 sorting of KNN values failure");
		}
		if(b>=m_fit.size() || b>=m_neighbours.size()){
			pagmo_throw(value_error,"SPEA2 sorting of KNN values failure");
		}
		double delta_a, delta_b;
		unsigned int i = 0;
		do{
			if(m_neighbours[b][i]>=m_fit.size() || m_neighbours[a][i]>=m_fit.size()){
				pagmo_throw(value_error,"SPEA2 sorting of KNN values failure");
			}
			delta_a = pagmo::util::neighbourhood::euclidian::distance(m_fit[a], m_fit[m_neighbours[a][i]]);
			delta_b = pagmo::util::neighbourhood::euclidian::distance(m_fit[b], m_fit[m_neighbours[b][i]]);
			i++;
		} while (i<m_neighbours[0].size() && delta_a == delta_b);
		return delta_a > delta_b;
	}
private:
	const std::vector<std::vector<pagmo::population::size_type>  > &m_neighbours;
	const std::vector<fitness_vector> &m_fit;
};

/// "Strength Pareto Evolutionary Algorithm (SPEA2)"
/**
 *
 * Strength Pareto Evolutionary Algorithm (SPEA2) is a multi-objective optimization algorithm.
 * The quality of an individual is measured taking into consideration its pareto strenght and its distance to its
 * K-th neighbour, where \f$ K=\sqrt{pop size + archive size} \f$
 * It uses an external archive in which are stored the non dominated solutions found so far.
 * The size of the archive is kept constant throughout the run by mean of a truncation operator taking into
 * consideration the distance of each individual to its closest neighbours.
 *
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 * @author Annalisa Riccardi (nina1983@gmail.com)
 *
 * @see Eckart Zitzler, Marco Laumanns, and Lothar Thiele -- "SPEA2: Improving the Strength Pareto Evolutionary Algorithm"
 **/
class __PAGMO_VISIBLE spea2: public base
{
public:
	spea2(int gen=100, double cr = 0.95, double eta_c = 10, double m = 0.01, double eta_m = 50, int archive_size = 0);
	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;

protected:
	std::string human_readable_extra() const;

private:
	struct spea2_individual {
		decision_vector x;
		fitness_vector f;
		constraint_vector c;
	};
	void compute_spea2_fitness(std::vector<double> &,
				int K,
				const std::vector<spea2_individual> &pop,
				const pagmo::problem::base &prob) const;
	std::vector<std::vector<population::size_type> > compute_domination_list(const pagmo::problem::base &,
																			const std::vector<fitness_vector> &,
																			const std::vector<constraint_vector> &) const;
	std::vector<population::size_type> compute_pareto_rank(const std::vector<std::vector<population::size_type> > &) const;
	pagmo::population::size_type tournament_selection(pagmo::population::size_type, pagmo::population::size_type,
													  const std::vector<population::size_type> &) const;
	std::vector<population::size_type> compute_domination_count(const std::vector<std::vector<population::size_type> > &) const;
	void crossover(decision_vector&, decision_vector&, pagmo::population::size_type, pagmo::population::size_type,
				   const std::vector<spea2_individual> &, const pagmo::problem::base &) const;
	void mutate(decision_vector&, const pagmo::problem::base&) const;
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
		ar & const_cast<int &>(m_archive_size);
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
	// Size of the archive
	int m_archive_size;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::spea2)

#endif // PAGMO_ALGORITHM_SPEA2_H
