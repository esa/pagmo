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

#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../types.h"
#include "base.h"
#include "mbh.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Specify an mbh algorithm with uniform neighbourhoods
 *
 * @param[in] local pagmo::algorithm to use as 'local' optimization method
 * @param[in] stop number of consecutive step allowed without any improvement
 * @param[in] perturb At the end of one iteration of mbh, each chromosome of each individual
 * will be perturbed within +-perturb, the same for the velocity. The integer part is treated the same way.
 * @throws value_error if stop is negative or perturb is negative
 */
mbh::mbh(const algorithm::base & local, int stop, double perturb):base(),m_stop(stop),m_perturb(1,perturb),m_screen_out(false)
{
	m_local = local.clone();
	if (stop < 0) {
		pagmo_throw(value_error,"number of consecutive step allowed without any improvement needs to be positive");
	}
	if (perturb < 0) {
		pagmo_throw(value_error,"perturb must be positive");
	}
}

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] local pagmo::algorithm to use as 'local' optimization method
 * @param[in] stop number of consecutive step allowed without any improvement
 * @param[in] perturb At the end of one iteration of mbh, the i-th chromosome of each individual
 * will be perturbed within +-perturb[i], the same for the velocity. The integer part is treated the same way.
 * @throws value_error if stop is negative or perturb[i] is negative
 */
mbh::mbh(const algorithm::base & local, int stop, const std::vector<double> &perturb):base(),m_stop(stop),m_perturb(perturb),m_screen_out(false)
{
	m_local = local.clone();
	if (stop < 0) {
		pagmo_throw(value_error,"number of consecutive step allowed without any improvement needs to be positive");
	}
	for (size_t i=0;i<perturb.size();++i)
	{
		if (perturb[i] < 0 ) {
			pagmo_throw(value_error,"perturb[.] must be positive");
		}
	}
	if (perturb.size()==0) pagmo_throw(value_error,"perturbation vector appears empty!!");
}

/// Copy constructor.
mbh::mbh(const mbh &algo):base(algo),m_local(algo.m_local->clone()),m_stop(algo.m_stop),m_perturb(algo.m_perturb),m_screen_out(algo.m_screen_out)
{}

/// Clone method.
base_ptr mbh::clone() const
{
	return base_ptr(new mbh(*this));
}

/// Evolve implementation.
/**
 * Run the MBH algorithm
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void mbh::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type D = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = pop.size();
	const problem::base::size_type Dc = D - prob_i_dimension;

	//Check if the perturbation vector has size 1, in which case it fills up the whole vector with
	//the same number
	if (m_perturb.size()==1)
	{
		for (problem::base::size_type i=1;i<D;++i) m_perturb.push_back(m_perturb[0]);
	}

	//Check that the perturbation vector size equals the size of the problem
	if (m_perturb.size()!=D)
	{
		pagmo_throw(value_error,"perturbation vector size does not match the problem size");
	}

	// Get out if there is nothing to do.
	if (m_stop == 0 || NP == 0) {
		return;
	}

	// Some dummies and temporary variables
	decision_vector tmp_x(D), tmp_v(D);
	double dummy, width;

	// Init the best fitness and constraint vector
	fitness_vector best_f = pop.get_individual(pop.get_best_idx()).cur_f;
	constraint_vector best_c = pop.get_individual(pop.get_best_idx()).cur_c;
	population best_pop(pop);

	int i = 0;

	//mbh main loop

	while (i<m_stop){

		//1. Perturb the current population (this could be moved in a pagmo::population method should other algorithm use it....
		for (population::size_type j =0; j < NP; ++j)
		{
			for (decision_vector::size_type k=0; k < Dc; ++k)
			{
				dummy = best_pop.get_individual(j).cur_x[k];
				width = m_perturb[k];
				tmp_x[k] = boost::uniform_real<double>(std::max(dummy-width,lb[k]),std::min(dummy+width,ub[k]))(m_drng);
				dummy = best_pop.get_individual(j).cur_v[k];
				tmp_v[k] = boost::uniform_real<double>(dummy-width,dummy+width)(m_drng);
			}

			for (decision_vector::size_type k=Dc; k < D; ++k)
			{
				dummy = best_pop.get_individual(j).cur_x[k];
				width = m_perturb[k];
				tmp_x[k] = boost::uniform_int<int>(std::max(dummy-width,lb[k]),std::min(dummy+width,ub[k]))(m_urng);
				dummy = best_pop.get_individual(j).cur_v[k];
				tmp_v[k] = boost::uniform_int<int>(std::max(dummy-width,lb[k]),std::min(dummy+width,ub[k]))(m_urng);
			}
			pop.set_x(j,tmp_x);
			pop.set_v(j,tmp_v);
		}

		//2. Evolve population with selected algorithm
		m_local->evolve(pop); i++;
		if (m_screen_out)
		{
			std::cout << i << ". " << "\tLocal solution: " << pop.get_individual(pop.get_best_idx()).cur_f << "\tGlobal best: " << best_f << std::endl;
		}

		//3. Reset counter if improved
		if (pop.problem().compare_fc(pop.get_individual(pop.get_best_idx()).cur_f,pop.get_individual(pop.get_best_idx()).cur_c,best_f,best_c) )
		{
			i = 0;
			best_f = pop.get_individual(pop.get_best_idx()).cur_f;
			best_c = pop.get_individual(pop.get_best_idx()).cur_c;
			if (m_screen_out) {
				std::cout << "New solution accepted. Constraints vector: " << best_c << '\n';
			}
			//update best population
			for (population::size_type j=0; j<pop.size();++j)
			{
				best_pop.set_x(j,pop.get_individual(j).cur_x);
				best_pop.set_v(j,pop.get_individual(j).cur_v);
			}
		}


	}
	//on exit set the population to the best one (discard perturbations)
	for (population::size_type j=0; j<pop.size();++j)
	{
		pop.set_x(j,best_pop.get_individual(j).cur_x);
		pop.set_v(j,best_pop.get_individual(j).cur_v);
	}
}

/// Activate screen output
/**
 * Activate screen output. Everytime a new champion is found the following information is printed
 * on the screen: Number of consecutive non improving iterations, best fitness
 *
 * @param[in] p true or false
 */
void mbh::screen_output(const bool p) {m_screen_out = p;}


/// Algorithm name
std::string mbh::get_name() const
{
	return "Generalized Monotonic Basin Hopping";
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string mbh::human_readable_extra() const
{
	std::ostringstream s;
	s << "algorithm:" << m_local->get_name() << ' ';
	s << "stop:" << m_stop << ' ';
	s << "perturb:" << m_perturb << ' ';
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::mbh);
