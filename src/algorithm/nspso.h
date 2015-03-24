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

#ifndef PAGMO_ALGORITHM_NSPSO_H
#define PAGMO_ALGORITHM_NSPSO_H

#include "../config.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace algorithm {

/// Non-dominated Sorting Particle Swarm Optimizer (NSPSO)
/**
 *
 * Non-dominated Sorting Particle Swarm Optimizer (NSPSO) is a modified version of PSO for multi-objective optimization.
 * It extends the basic ideas of PSO by making a better use of personal bests and offspring for non-dominated comparison.
 * In order to increase the diversity of the pareto front it is possible to choose between 3 different niching methods:
 * crowding distance, niche count and maxmin.
 *
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 * @author Annalisa Riccardi (nina1983@gmail.com)
 *
 * @see "Xiaodong Li - A Non-dominated Sorting Particle Swarm Optimizer for Multiobjective Optimization"
 * @see "Xiaodong Li - Better Spread and Convergence: Particle Swarm Multiobjective Optimization Using the Maximin Fitness Function"
 * @see "Carlos M. Fonseca, Peter J. Fleming - Genetic Algorithms for Multiobjective Optimization: Formulation, Discussion and Generalization"
 **/


class __PAGMO_VISIBLE nspso: public base
{
public:
	/// Mechanism used to asses diversity
	enum diversity_mechanism_type {
	   CROWDING_DISTANCE=0, ///< The crowding distance is used
	   NICHE_COUNT=1, ///< The ....
	   MAXMIN=2 ///< The ....
	};
	nspso(int gen=100, double minW = 0.4, double maxW = 1.0, double C1 = 2.0, double C2 = 2.0,
		  double CHI = 1.0, double v_coeff = 0.5, int leader_selection_range = 10, diversity_mechanism_type = CROWDING_DISTANCE);

	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;

protected:
	std::string human_readable_extra() const;

private:
	struct one_dim_fit_comp {
		one_dim_fit_comp(const std::vector<fitness_vector> &fit, fitness_vector::size_type dim):m_fit(fit),m_dim(dim){};
		bool operator()(const population::size_type& idx1, const population::size_type& idx2) const
		{
			return m_fit[idx1][m_dim] < m_fit[idx2][m_dim];
		}
		const std::vector<fitness_vector>& m_fit;
		fitness_vector::size_type m_dim;
	};

	struct crowding_pareto_comp{
		crowding_pareto_comp(const std::vector<population::size_type> &pareto_rank, const std::vector<double> &crowding_d):m_pareto_rank(pareto_rank),m_crowding_d(crowding_d){};
		bool operator()(const population::size_type& idx1, const population::size_type& idx2) const
		{
			if (m_pareto_rank[idx1] == m_pareto_rank[idx2]) {
				return (m_crowding_d[idx1] > m_crowding_d[idx2]);
			}
			else {
				return (m_pareto_rank[idx1] < m_pareto_rank[idx2]);
			}
		}
		const std::vector<population::size_type> &m_pareto_rank;
		const std::vector<double> &m_crowding_d;
	};

	struct nspso_individual {
		decision_vector cur_x;
		decision_vector best_x;
		decision_vector cur_v;
		fitness_vector cur_f;
		fitness_vector best_f;
		constraint_vector cur_c;
		constraint_vector best_c;
	};


	double minfit(unsigned int, unsigned int, const std::vector<fitness_vector> &) const;
	void compute_maxmin(std::vector<double> &, const std::vector<fitness_vector> &) const;
	void compute_niche_count(std::vector<int> &, const std::vector<std::vector<double> > &, double) const;
	double euclidian_distance(const std::vector<double> &, const std::vector<double> &) const;
	std::vector<std::vector<population::size_type> > compute_domination_list(const pagmo::problem::base &,
																			const std::vector<fitness_vector> &,
																			const std::vector<constraint_vector> &) const;
	std::vector<population::size_type> compute_domination_count(const std::vector<std::vector<population::size_type> > &) const;
	std::vector<population::size_type> compute_pareto_rank(const std::vector<std::vector<population::size_type> > &) const;
	std::vector<std::vector<population::size_type> > compute_pareto_fronts(const std::vector<population::size_type> &) const;
	std::vector<std::vector<population::size_type> > compute_pareto_fronts(const pagmo::problem::base &,
																		   const std::vector<fitness_vector> &,
																		   const std::vector<constraint_vector> &) const;
	std::vector<double> compute_crowding_d(const std::vector<fitness_vector> &, const std::vector<std::vector<population::size_type> > &) const;
	fitness_vector compute_ideal(const std::vector<fitness_vector> &, const std::vector<population::size_type> &) const;
	fitness_vector compute_nadir(const std::vector<fitness_vector> &, const std::vector<population::size_type> &) const;

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<int &>(m_gen);
		ar & const_cast<double &>(m_minW);
		ar & const_cast<double &>(m_maxW);
		ar & const_cast<double &>(m_C1);
		ar & const_cast<double &>(m_C2);
		ar & const_cast<double &>(m_CHI);
		ar & const_cast<double &>(m_v_coeff);
		ar & const_cast<int &>(m_leader_selection_range);
		ar & const_cast<diversity_mechanism_type &>(m_diversity_mechanism);

	}
	//Number of generations
	const int m_gen;
	const double m_minW;
	const double m_maxW;
	const double m_C1;
	const double m_C2;
	const double m_CHI;
	const double m_v_coeff;
	const int m_leader_selection_range;
	const diversity_mechanism_type m_diversity_mechanism;

};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::nspso)

#endif // PAGMO_ALGORITHM_NSPSO_H
