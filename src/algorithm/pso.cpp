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

#include "pso.h"
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>


namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] gen number of generations.
 * @param[in] m_omega particle inertia
 * @param[in] eta1 cognitive component of the particle
 * @param[in] eta2 social component of the particle
 * @param[in] vcoeff velocity coefficient (determining the maximum allowed particle velocity)
 * @param[in] variant algorithm variant to use
 * @throws value_error if m_omega is not in the [0,1] interval, eta1, eta2 are not in the [0,1] interval,
 * vcoeff is not in ]0,1[, variant is not one of 1 .. 4
 */
pso::pso(int gen, double m_omega, double eta1, double eta2, double vcoeff, int variant):base(),m_gen(gen),m_omega(m_omega),m_eta1(eta1),m_eta2(eta2),m_vcoeff(vcoeff),m_variant(variant) {
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}

	if (m_omega < 0 || m_omega > 1) {
		pagmo_throw(value_error,"the particles inertia must be in the [0,1] range");
	}

	if (eta1 < 0 || eta2 < 0 || eta1 > 4 || eta2 > 4) {
		pagmo_throw(value_error,"the eta parameters must be in the [0,1] range");
	}

	if (vcoeff <= 0 || vcoeff >= 1) {
		pagmo_throw(value_error,"algorithm variant must be one of 1 ... 4");
	}

	if (variant < 1 || variant > 4) {
		pagmo_throw(value_error,"algorithm variant must be one of 1 ... 4");
	}


}

/// Clone method.
base_ptr pso::clone() const
{
	return base_ptr(new pso(*this));
}

/// Evolve implementation.
/**
 * Run the PSO algorithm for the number of generations specified in the constructors.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void pso::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type D = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension(), prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = pop.size();
	const problem::base::size_type Dc = D - prob_i_dimension;


	//We perform some checks to determine wether the problem/population are suitable for PSO
	if ( Dc == 0 ) {
		pagmo_throw(value_error,"There is no continuous part in the problem decision vector for PSO to optimise");
	}

	if ( prob_c_dimension != 0 ) {
		pagmo_throw(value_error,"The problem is not box constrained and PSO is not suitable to solve it");
	}

	if ( prob_f_dimension != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and PSO is not suitable to solve it");
	}

	// Get out if there is nothing to do.
	if (NP == 0 || m_gen == 0) {
		return;
	}

	// Some vectors used during evolution are allocated here.
	std::vector<double> dummy(D,0);			//used for initialisation purposes
	std::vector<decision_vector > X(NP,dummy);	//particle position
	std::vector<decision_vector > V(NP,dummy);	//particle velocity
	std::vector<fitness_vector> fit(NP);		//particle fitness

	fitness_vector gbfit;				//global best fitness
	decision_vector gbX(D);				//global best chromosome

	std::vector<fitness_vector> lbfit(NP);		//personal best fitness
	std::vector<decision_vector> lbX(NP,dummy);	//personal best chromosome
	decision_vector minv(D),maxv(D);		//Maximum and minumum velocity allowed

	double vwidth;					//Temporary variable


	// Copy the particle positions, their velocities and their fitness
	for ( population::size_type i = 0; i<NP; i++ ) {
		X[i]	=	pop.get_individual(i).cur_x;
		V[i]	=	pop.get_individual(i).cur_v;
		fit[i]	=	pop.get_individual(i).cur_f;
	}

	// Initialise the minimum and maximum velocity
	for ( problem::base::size_type i = 0; i<D; i++ ) {
		vwidth = (ub[i]-lb[i]) * m_vcoeff;
		minv[i] = -1.0*vwidth;
		maxv[i] = vwidth;
	}

	// Initialise the global and local bests
	gbX=pop.champion().x;
	gbfit=pop.champion().f;

	for (population::size_type i=0; i < NP;++i){
	lbX[i] = pop.get_individual(i).best_x;
	lbfit[i] = pop.get_individual(i).best_f;
	}

	double r1,r2 = 0;
	// Main PSO loop
	for (int j = 0; j < m_gen; ++j) {
		//For each particle in the swarm
		for (population::size_type ii = 0; ii< NP; ii++) {

			/*-------PSO canonical--------------------------------------------------------------------*/
			/*-------The classical PSO velocity update startegy---------------------------------------*/
			if (m_variant==1) {
				r1 = m_drng();
				r2 = m_drng();
				for (problem::base::size_type jj = 0; jj< Dc; jj++) {
					V[ii][jj] = m_omega * V[ii][jj] + m_eta1 * r1 * (lbX[ii][jj] - X[ii][jj]) + m_eta2 * r2 * (gbX[jj] - X[ii][jj]);
				}
			}

			/*-------PSO canonical with equal random weights of social and cognitive components-------*/
			/*-------In our experience few problems benefit a lot from having r1=r2. You may check----*/
			/*-------with Rastrigin-------------------------------------------------------------------*/
			if (m_variant==2) {
				r1 = m_drng();
				for (problem::base::size_type jj = 0; jj< Dc; jj++) {
					V[ii][jj] = m_omega * V[ii][jj] + m_eta1 * r1 * (lbX[ii][jj] - X[ii][jj]) + m_eta2 * r1 * (gbX[jj] - X[ii][jj]);
				}
			}

			/*-------Variant of PSO strategy 1 that has r1 and r2 randomly generated for each----------*/
			/*-------component. This is also the version that was implemented--------------------------*/
			/*-------in early versions of pagmo--------------------------------------------------------*/
			if (m_variant==3) {
				for (problem::base::size_type jj = 0; jj< Dc; jj++) {
					r1 = m_drng();
					r2 = m_drng();
					V[ii][jj] = m_omega * V[ii][jj] + m_eta1 * r1 * (lbX[ii][jj] - X[ii][jj]) + m_eta2 * r2 * (gbX[jj] - X[ii][jj]);
				}
			}

			/*-------Variant of PSO strategy 2 that has r1 randomly generated for each-----------------*/
			/*-------component.------------------------------------------------------------------------*/
			if (m_variant==4) {
				for (problem::base::size_type jj = 0; jj< Dc; jj++) {
					r1 = m_drng();
					V[ii][jj] = m_omega * V[ii][jj] + m_eta1 * r1 * (lbX[ii][jj] - X[ii][jj]) + m_eta2 * r1 * (gbX[jj] - X[ii][jj]);
				}
			}

			//We now check that the velocity does not exceed the maximum allowed per component
			//and we perform the position update and the feasibility correction
			for (problem::base::size_type jj = 0; jj< Dc; jj++) {

				if ( V[ii][jj] > maxv[jj] )
					V[ii][jj] = maxv[jj];

				else if ( V[ii][jj] < minv[jj] )
					V[ii][jj] = minv[jj];

				//update position
				X[ii][jj] = X[ii][jj] + V[ii][jj];

				//feasibility correction
				if (X[ii][jj] < lb[jj])
					X[ii][jj] = boost::uniform_real<double>(lb[jj],ub[jj])(m_drng);

				else if (X[ii][jj] > ub[jj])
					X[ii][jj] = boost::uniform_real<double>(lb[jj],ub[jj])(m_drng);
			}

			//We evaluate here the new individual fitness as to be able to update the global best in real time
			prob.objfun(fit[ii],X[ii]);
			//update local and global best
			if ( prob.compare_fitness(fit[ii],lbfit[ii])) {
				lbfit[ii] = fit[ii];	//local best
				lbX[ii] = X[ii];
				pop.set_x(ii,X[ii]);
				pop.set_v(ii,V[ii]);
				if ( prob.compare_fitness(fit[ii],gbfit) ) {
					gbfit = fit[ii];	//global best
					gbX = X[ii];
				}
			}
		} //End of loop on the population members
	} // end of main PSO loop
}

/// Algorithm name
std::string pso::get_name() const
{
	return "Particle Swarm optimization";
}


/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string pso::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	s << "omega:" << m_omega << ' ';
	s << "eta1:" << m_eta1 << ' ';
	s << "eta2:" << m_eta2 << ' ';
	s << "variant:" << m_variant << ' ';
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::pso);
