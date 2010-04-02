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
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] local pagmo::algorithm to use as 'local' optimization method
 * @param[in] stop number of consecutive step allowed without any improvement
 * @param[in] perturb At the end of one iteration of mbh, each chromosome of each population individual
 * will be perturbed by +-perturb*(ub - lb)
 * @throws value_error if stop is negative or perturb is not in ]0,1]
 */
mbh::mbh(base local, int, stop, double perturb):base(),m_stop(stop),m_perturb(perturb)
{
	m_local = local.clone();
	if (stop < 0) {
		pagmo_throw(value_error,"number of consecutive step allowed without any improvement needs to be positive");
	}
	if (perturb <= 0 || perturb > 1) {
		pagmo_throw(value_error,"perturb must be in ]0,1]");
	}
}

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
	const problem::base::size_type D = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension(), prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = pop.size();
	const problem::base::size_type Dc = D - prob_i_dimension;

	// Get out if there is nothing to do.
	if (stop == 0 || NP == 0) {
		return;
	}

	// Init the best chromosome, fitness
	decision_vector best = pop.champion().x;
	fitness_vector current = pop.champion().f;

	int i = 0;

	//mbh main loop

	while (i<stop){
		pop = m_local.evolve(pop); i++;
		if (prob.compare_fitness(pop.champion().f),current){
			i = 0;
			best = pop[0].getDecisionVector();
		}
		//perturb the population

  49                 for (size_t i=0;i<deme.problem().getDimension();i++){

  50                         current[i] = best[i] + (drng()*2 -1)*perturb*(ub[i]-lb[i]);

  51                         if (current[i] > ub[i]) current[i] = ub[i] - drng() * perturb*(ub[i]-lb[i]);

  52                         if (current[i] < lb[i]) current[i] = lb[i] + drng() * perturb*(ub[i]-lb[i]);;

  53                 }



}



/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string mbh::human_readable_extra() const
{
	std::ostringstream s;
	s << "\tGenerations:\t" << m_gen << '\n';
	s << "\tWeight parameter (F):\t\t" << m_f << '\n';
	s << "\tCrossover parameter (CR):\t" << m_cr << '\n';
	s << "\tStrategy selected:\t" << m_strategy << '\n';
	return s.str();
}

}} //namespaces
