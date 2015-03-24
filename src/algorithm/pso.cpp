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

#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <vector>
#include <cmath>
#include <iostream>

#include "pso.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] gen number of generations
 * @param[in] omega particles' inertia weight, or alternatively, the constriction coefficient (usage depends on the variant used)
 * @param[in] eta1 magnitude of the force, applied to the particle's velocity, in the direction of its previous best position
 * @param[in] eta2 magnitude of the force, applied to the particle's velocity, in the direction of the best position in its neighborhood
 * @param[in] vcoeff velocity coefficient (determining the maximum allowed particle velocity)
 * @param[in] variant algorithm variant to use
 * @param[in] neighb_type swarm topology to use
 * @param[in] neighb_param if the lbest topology is selected (neighb_type=2), it represents each particle's indegree
 * (also outdegree) in the swarm topology. Particles have neighbours up
 * to a radius of k = neighb_param / 2 in the ring. If the Randomly-varying neighbourhood topology
 * is selected (neighb_type=4), neighb_param represents each particle's maximum outdegree in the swarm topology.
 * The minimum outdegree is 1 (the particle always connects back to itself).
 * @throws value_error if m_omega is not in the [0,1] interval, eta1, eta2 are not in the [0,1] interval,
 * vcoeff is not in ]0,1], variant is not one of 1 .. 6, neighb_type is not one of 1 .. 4
 */
pso::pso(int gen, double omega, double eta1, double eta2, double vcoeff, int variant, int neighb_type, int neighb_param):base(),m_gen(gen),m_omega(omega),m_eta1(eta1),m_eta2(eta2),m_vcoeff(vcoeff),m_variant(variant),m_neighb_type(neighb_type),m_neighb_param(neighb_param) {
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	
	if (m_omega < 0 || m_omega > 1) {
		if( variant < 5 )
			// variants using Inertia weight
			pagmo_throw(value_error,"the particles' inertia must be in the [0,1] range");
		else
			// variants using Constriction coefficient
			pagmo_throw(value_error,"the particles' constriction coefficient must be in the [0,1] range");
	}
	
	if (eta1 < 0 || eta2 < 0 || eta1 > 4 || eta2 > 4) {
		pagmo_throw(value_error,"the eta parameters must be in the [0,4] range");
	}
	
	if (vcoeff <= 0 || vcoeff > 1) {
		pagmo_throw(value_error,"fraction of variables' range in which velocity may vary should be in ]0,1]");
	}
	
	if (variant < 1 || variant > 6) {
		pagmo_throw(value_error,"algorithm variant must be one of 1 ... 6");
	}
	
	if (neighb_type < 1 || neighb_type > 4) {
		pagmo_throw(value_error,"swarm topology variant must be one of 1 ... 4");
	}
	// And neighb param???
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
	const problem::base             &prob = pop.problem();
	const problem::base::size_type   D = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension(), prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();
	const problem::base::size_type   Dc = D - prob_i_dimension;
	const decision_vector           &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type      swarm_size = pop.size();
	
	
	//We perform some checks to determine wether the problem/population are suitable for PSO
	if( Dc == 0 ){
		pagmo_throw(value_error,"There is no continuous part in the problem decision vector for PSO to optimise");
	}
	
	if( prob_c_dimension != 0 ){
		pagmo_throw(value_error,"The problem is not box constrained and PSO is not suitable to solve it");
	}
	
	if( prob_f_dimension != 1 ){
		pagmo_throw(value_error,"The problem is not single objective and PSO is not suitable to solve it");
	}
	
	// Get out if there is nothing to do.
	if (swarm_size == 0 || m_gen == 0) {
		return;
	}

	// Some vectors used during evolution are allocated here.
	std::vector<double> dummy(Dc,0);				// used for initialisation purposes
	
	std::vector<decision_vector>  X(swarm_size,dummy);	// particles' current positions
	std::vector<fitness_vector>   fit(swarm_size);		// particles' current fitness values
	
	std::vector<decision_vector > V(swarm_size,dummy);	// particles' velocities
	
	std::vector<decision_vector> lbX(swarm_size,dummy);	// particles' previous best positions
	std::vector<fitness_vector>  lbfit(swarm_size);		// particles' fitness values at their previous best positions
	
	std::vector< std::vector<int> > neighb(swarm_size);	// swarm topology (iterators over indexes of each particle's neighbors in the swarm)
	
	decision_vector best_neighb(Dc);			// search space position of particles' best neighbor
	fitness_vector best_fit;				// fitness at the best found search space position (tracked only when using topologies 1 or 4)
	bool best_fit_improved;					// flag indicating whether the best solution's fitness improved (tracked only when using topologies 1 or 4)
	
	
	decision_vector minv(Dc), maxv(Dc);			// Maximum and minimum velocity allowed
	
	double vwidth;						// Temporary variable
	double new_x;						// Temporary variable
	
	population::size_type    p;		// for iterating over particles
	population::size_type    n;		// for iterating over particles's neighbours
	problem::base::size_type d;		// for iterating over problem dimensions
	
	
	// Initialise the minimum and maximum velocity
	for( d = 0; d < Dc; d++ ){
		vwidth  = ( ub[d] - lb[d] ) * m_vcoeff;
		minv[d] = -1.0 * vwidth;
		maxv[d] = vwidth;
	}
	
	// Copy the particle positions, their velocities and their fitness
	for( p = 0; p < swarm_size; p++ ){
		X[p]   = pop.get_individual(p).cur_x;
		V[p]   = pop.get_individual(p).cur_v;
		fit[p] = pop.get_individual(p).cur_f;
	}
	
	// Initialise particles' previous best positions
	for( p = 0; p < swarm_size; p++ ){
		lbX[p]   = pop.get_individual(p).best_x;
		lbfit[p] = pop.get_individual(p).best_f;
	}
	
	// Initialize the Swarm's topology
	switch( m_neighb_type ){
		case 1:  initialize_topology__gbest( pop, best_neighb, best_fit, neighb ); break;
		case 3:  initialize_topology__von( neighb ); break;
		case 4:  initialize_topology__adaptive_random( neighb );
			best_fit = pop.champion().f;	// need to track improvements in best found fitness, to know when to rewire
			break;
		case 2:
		default: initialize_topology__lbest( neighb );
	}
	
	
	// auxiliary varibables specific to the Fully Informed Particle Swarm variant
	double acceleration_coefficient = m_eta1 + m_eta2;
	double sum_forces;
	
	double r1 = 0.0;
	double r2 = 0.0;
	
	/* --- Main PSO loop ---
	 */
	// For each generation
	for( int g = 0; g < m_gen; ++g ){
		
		best_fit_improved = false;
		
		// For each particle in the swarm
		for( p = 0; p < swarm_size; p++ ){
			
			// identify the current particle's best neighbour
			// . not needed if m_neighb_type == 1 (gbest): best_neighb directly tracked in this function
			// . not needed if m_variant == 6 (FIPS): all neighbours are considered, no need to identify the best one
			if( m_neighb_type != 1 && m_variant != 6)
				best_neighb = particle__get_best_neighbor( p, neighb, lbX, lbfit, prob );
				
			
			/*-------PSO canonical (with inertia weight) ---------------------------------------------*/
			/*-------Original algorithm used in PaGMO paper-------------------------------------------*/
			if( m_variant == 1 ){
				for( d = 0; d < Dc; d++ ){
					r1 = m_drng();
					r2 = m_drng();
					V[p][d] = m_omega * V[p][d] + m_eta1 * r1 * (lbX[p][d] - X[p][d]) + m_eta2 * r2 * (best_neighb[d] - X[p][d]);
				}
			}
			
			/*-------PSO canonical (with inertia weight) ---------------------------------------------*/
			/*-------and with equal random weights of social and cognitive components-----------------*/
			/*-------Check with Rastrigin-------------------------------------------------------------*/
			else if( m_variant == 2 ){
				for( d = 0; d < Dc; d++ ){
					r1 = m_drng();
					V[p][d] = m_omega * V[p][d] + m_eta1 * r1 * (lbX[p][d] - X[p][d]) + m_eta2 * r1 * (best_neighb[d] - X[p][d]);
				}
			}
			
			/*-------PSO variant (commonly mistaken in literature for the canonical)----------------*/
			/*-------Same random number for all components------------------------------------------*/
			else if( m_variant == 3 ){
				r1 = m_drng();
				r2 = m_drng();
				for( d = 0; d < Dc; d++ ){
					V[p][d] = m_omega * V[p][d] + m_eta1 * r1 * (lbX[p][d] - X[p][d]) + m_eta2 * r2 * (best_neighb[d] - X[p][d]);
				}
			}
			
			/*-------PSO variant (commonly mistaken in literature for the canonical)----------------*/
			/*-------Same random number for all components------------------------------------------*/
			/*-------and with equal random weights of social and cognitive components---------------*/
			else if( m_variant == 4 ){
				r1 = m_drng();
				for( d = 0; d < Dc; d++ ){
					V[p][d] = m_omega * V[p][d] + m_eta1 * r1 * (lbX[p][d] - X[p][d]) + m_eta2 * r1 * (best_neighb[d] - X[p][d]);
				}
			}
			
			/*-------PSO variant with constriction coefficients------------------------------------*/
			/*  ''Clerc's analysis of the iterative system led him to propose a strategy for the
			 *  placement of "constriction coefficients" on the terms of the formulas; these
			 *  coefficients controlled the convergence of the particle and allowed an elegant and
			 *  well-explained method for preventing explosion, ensuring convergence, and
			 *  eliminating the arbitrary Vmax parameter. The analysis also takes the guesswork
			 *  out of setting the values of phi_1 and phi_2.''
			 *  ''this is the canonical particle swarm algorithm of today.''
			 *  [Poli et al., 2007] http://dx.doi.org/10.1007/s11721-007-0002-0
			 *  [Clerc and Kennedy, 2002] http://dx.doi.org/10.1109/4235.985692
			 *  
			 *  This being the canonical PSO of today, this variant is set as the default in PaGMO.
			 *-------------------------------------------------------------------------------------*/
			else if( m_variant == 5 ){
				for( d = 0; d < Dc; d++ ){
					r1 = m_drng();
					r2 = m_drng();
					V[p][d] = m_omega * ( V[p][d] + m_eta1 * r1 * (lbX[p][d] - X[p][d]) + m_eta2 * r2 * (best_neighb[d] - X[p][d]) );
				}
			}
			
			/*-------Fully Informed Particle Swarm-------------------------------------------------*/
			/*  ''Whereas in the traditional algorithm each particle is affected by its own
			 *  previous performance and the single best success found in its neighborhood, in
			 *  Mendes' fully informed particle swarm (FIPS), the particle is affected by all its
			 *  neighbors, sometimes with no influence from its own previous success.'' 
			 *  ''With good parameters, FIPS appears to find better solutions in fewer iterations
			 *  than the canonical algorithm, but it is much more dependent on the population topology.''
			 *  [Poli et al., 2007] http://dx.doi.org/10.1007/s11721-007-0002-0
			 *  [Mendes et al., 2004] http://dx.doi.org/10.1109/TEVC.2004.826074
			 *-------------------------------------------------------------------------------------*/
			else if( m_variant == 6 ){
				for( d = 0; d < Dc; d++ ){
					sum_forces = 0.0;
					for( n = 0; n < neighb[p].size(); n++ )
						sum_forces += m_drng() * acceleration_coefficient * ( lbX[ neighb[p][n] ][d] - X[p][d] );
					
					V[p][d] = m_omega * ( V[p][d] + sum_forces / neighb[p].size() );
				}
			}
			
			// We now check that the velocity does not exceed the maximum allowed per component
			// and we perform the position update and the feasibility correction
			for( d = 0; d < Dc; d++ ){
				
				if( V[p][d] > maxv[d] )
					V[p][d] = maxv[d];
				
				else if( V[p][d] < minv[d] )
					V[p][d] = minv[d];
				
				// update position
				new_x = X[p][d] + V[p][d];
				
				// feasibility correction
				// (velocity updated to that which would have taken the previous position
				// to the newly corrected feasible position)
				if( new_x < lb[d] ){
					new_x = lb[d];
					V[p][d] = 0.0;
//					new_x = boost::uniform_real<double>(lb[d],ub[d])(m_drng);
//					V[p][d] = new_x - X[p][d];
				}
				else if( new_x > ub[d] ){
					new_x = ub[d];
					V[p][d] = 0.0;
//					new_x = boost::uniform_real<double>(lb[d],ub[d])(m_drng);
//					V[p][d] = new_x - X[p][d];
				}
				
				X[p][d] = new_x;
			}
			
			// We evaluate here the new individual fitness as to be able to update the global best in real time
			prob.objfun( fit[p], X[p] );
			m_fevals++;
			
			if( prob.compare_fitness( fit[p], lbfit[p] ) ){
				// update the particle's previous best position
				lbfit[p] = fit[p];
				lbX[p] = X[p];
				
				// update the best position observed so far by any particle in the swarm
				// (only performed if swarm topology is gbest)
				if( ( m_neighb_type == 1 || m_neighb_type == 4 ) && prob.compare_fitness( fit[p], best_fit ) ){
					best_neighb = X[p];
					best_fit    = fit[p];
					best_fit_improved = true;
				}
			}
		
		} //End of loop on the population members
		
		// reset swarm topology if no improvement was observed in the best found fitness value
		if( m_neighb_type == 4 && !best_fit_improved )
			initialize_topology__adaptive_random( neighb );
		
	} // end of main PSO loop
	
	
	// copy particles' positions & velocities back to the main population
	for( p = 0; p < swarm_size; p++ ){
		pop.set_x( p, lbX[p] );		// sets: cur_x, cur_f, best_x, best_f
		pop.set_x( p, X[p] );		// sets: cur_x, cur_f
		pop.set_v( p, V[p] );		// sets: cur_v
	}
}


/**
 *  @brief Get information on the best position already visited by any of a particle's neighbours
 *  
 *  @param[in] pidx index to the particle under consideration
 *  @param[in] neighb definition of the swarm's topology
 *  @param[in] lbX particles' previous best positions
 *  @param[in] lbfit particles' fitness values at their previous best positions
 *  @param[in] prob problem undergoing optimization
 *  @return best position already visited by any of the considered particle's neighbours
 */
decision_vector pso::particle__get_best_neighbor( population::size_type pidx, std::vector< std::vector<int> > &neighb, const std::vector<decision_vector> &lbX, const std::vector<fitness_vector> &lbfit, const problem::base &prob ) const
{
	population::size_type  nidx, bnidx;		// neighbour index; best neighbour index
	
	switch( m_neighb_type ){
		case 1: // { gbest }
			// ERROR: execution should not reach this point, as the global best position is not tracked using the neighb vector
			pagmo_throw(value_error,"particle__get_best_neighbor() invoked while using a gbest swarm topology");
			break;
		case 2: // { lbest }
		case 3: // { von }
		case 4: // { adaptive random }
		default:
			// iterate over indexes of the particle's neighbours, and identify the best
			bnidx = neighb[pidx][0];
			for( nidx = 1; nidx < neighb[pidx].size(); nidx++ )
				if( prob.compare_fitness( lbfit[ neighb[pidx][nidx] ], lbfit[ bnidx ] ) )
					bnidx = neighb[pidx][nidx];
			return lbX[bnidx];
	}
}


/**
 *  @brief Defines the Swarm topology as a fully connected graph, where particles are influenced by all other particles in the swarm
 *  
 *  ''The earliest reported particle swarm version [3], [4] used a kind of
 *  topology that is known as gbest. The source of social influence on each
 *  particle was the best performing individual in the entire population. This
 *  is equivalent to a sociogram or social network where every individual is
 *  connected to every other one.'' \n
 *  ''The gbest topology (i.e., the biggest neighborhood possible) has often
 *  been shown to converge on optima more quickly than lbest, but gbest is also
 *  more susceptible to the attraction of local optima since the population
 *  rushes unanimously toward the first good solution found.'' \n
 *  [Kennedy and Mendes, 2006] http://dx.doi.org/10.1109/TSMCC.2006.875410
 *  
 *  @param[in] pop pagmo::population being evolved
 *  @param[out] gbX best search space position already visited by the swarm
 *  @param[out] gbfit best fitness value in the swarm
 *  @param[out] neighb definition of the swarm's topology
 */
void pso::initialize_topology__gbest( const population &pop, decision_vector &gbX, fitness_vector &gbfit, std::vector< std::vector<int> > &neighb ) const
{
	// The best position already visited by the swarm will be tracked in pso::evolve() as particles are evaluated.
	// Here we define the initial values of the variables that will do that tracking.
	gbX   = pop.champion().x;
	gbfit = pop.champion().f;
	
	/* The usage of a gbest swarm topology along with a FIPS (fully informed particle swarm) velocity update formula
	 * is discouraged. However, because a user might still configure such a setup, we must ensure FIPS has access to
	 * the list of indices of particles' neighbours:
	 */
	if( m_variant == 6 ){
		unsigned int i;
		for( i = 0; i < neighb.size(); i++ )
			neighb[0].push_back( i );
		for( i = 1; i < neighb.size(); i++ )
			neighb[i] = neighb[0];
	}
}


/**
 *  @brief Defines the Swarm topology as a ring, where particles are influenced by their immediate neighbors to either side
 *  
 *  ''The L-best-k topology consists of n nodes arranged in a ring, in which
 *  node i is connected to each node in {(i+j) mod n : j = +-1,+-2, . . . ,+-k}.'' \n
 *  [Mohais et al., 2005] http://dx.doi.org/10.1007/11589990_80
 *  
 *  neighb_param represents each particle's indegree (also outdegree) in the swarm topology.
 *	Particles have neighbours up to a radius of k = neighb_param / 2 in the ring.
 *  
 *  @param[out] neighb definition of the swarm's topology
 */
void pso::initialize_topology__lbest( std::vector< std::vector<int> > &neighb ) const
{
	int swarm_size = neighb.size();
	int pidx;		// for iterating over particles
	int nidx, j;	// for iterating over particles' neighbours
	
	int radius = m_neighb_param / 2;
	
	for( pidx = 0; pidx < swarm_size; pidx++ ){
		for( j = -radius; j <= radius; j++ ){
			if( j == 0 ) j++;
			nidx = (pidx + j) % swarm_size;
			if( nidx < 0 ) nidx = swarm_size + nidx;
			neighb[pidx].push_back( nidx );
		}
	}
}


/*! @brief Von Neumann neighborhood
 *  (increments on particles' lattice coordinates that produce the coordinates of their neighbors)
 *  
 *  The von Neumann neighbourhood of a point includes all the points at a Hamming distance of 1.
 *  
 *  - http://en.wikipedia.org/wiki/Von_Neumann_neighborhood
 *  - http://mathworld.wolfram.com/vonNeumannNeighborhood.html
 *  - http://en.wikibooks.org/wiki/Cellular_Automata/Neighborhood
 */
const int	vonNeumann_neighb_diff[4][2] = { {-1,0}, {1,0}, {0,-1}, {0,1} };

/**
 *  @brief Arranges particles in a lattice, where each interacts with its immediate 4 neighbors to the N, S, E and W.
 *  
 *  ''The population is arranged in a rectangular matrix, for instance, 5 x 4
 *  in a population of 20 individuals, and each individual is connected to
 *  the individuals above, below and on both of its sides, wrapping the edges'' \n
 *  [Kennedy and Mendes, 2006] http://dx.doi.org/10.1109/TSMCC.2006.875410
 *  
 *  ''Given a population of size n, the von Neumann neighborhood was
 *  configured into r rows and c columns, where r is the smallest integer
 *  less than or equal to sqrt(n) that evenly divides n and c = n / r'' \n
 *  [Mohais et al., 2005] http://dx.doi.org/10.1007/11589990_80 \n
 *  (there's an error in the description above: "smallest integer" should
 *  instead be "highest integer")
 *  
 *  @param[out] neighb definition of the swarm's topology
 */
void pso::initialize_topology__von( std::vector< std::vector<int> > &neighb ) const
{
	int swarm_size = neighb.size();
	int	cols, rows;		// lattice structure
	int	pidx, nidx;		// particle and neighbour indices, in the swarm and neighbourhood vectors
	int	p_x, p_y;		// particle's coordinates in the lattice
	int	n_x, n_y;		// particle neighbor's coordinates in the lattice
	
	rows = std::sqrt( swarm_size );
	while( swarm_size % rows != 0 )
		rows -= 1;
	cols = swarm_size / rows;
	
	for( pidx = 0; pidx < swarm_size; pidx++ ){
		p_x = pidx % cols;
		p_y = pidx / cols;
		
		for( nidx = 0; nidx < 4; nidx++ ){
			n_x = ( p_x + vonNeumann_neighb_diff[nidx][0] ) % cols;  if( n_x < 0 ) n_x = cols + n_x;	// sign of remainder(%) in a division when at least one of the operands is negative is compiler implementation specific. The 'if' here ensures the same behaviour across compilers
			n_y = ( p_y + vonNeumann_neighb_diff[nidx][1] ) % rows;  if( n_y < 0 ) n_y = rows + n_y;
			
			neighb[pidx].push_back( n_y * cols + n_x );
		}
	}
}


/**
 *  @brief Arranges particles in a random graph having a parameterized maximum outdegree; the graph changes adaptively over time
 *  
 *	''At the very beginning, and after each unsuccessful iteration (no
 *	improvement of the best known fitness value), the graph of the information
 *	links is modified: each particle informs at random K particles (the same
 *	particle may be chosen several times), and informs itself. The parameter K
 *	is usually set to 3. It means that each particle informs at less one
 *	particle (itself), and at most K+1 particles (including itself). It also
 *	means that each particle can be informed by any number of particles between
 *	1 and S. However, the distribution of the possible number of "informants"
 *	is not uniform. On average, a particle is often informed by about K others,
 *	but it may be informed by a much larger number of particles with a small
 *	probability'' \n
 *  [Maurice Clerc, 2011] Standard Particle Swarm Optimisation, From 2006 to 2011 \n
 *	http://clerc.maurice.free.fr/pso/SPSO_descriptions.pdf
 *	
 *  neighb_param represents each particle's maximum outdegree in the swarm topology.
 *	The minimum outdegree is 1 (the particle always connects back to itself).
 *  
 *  @param[out] neighb definition of the swarm's topology
 */
void pso::initialize_topology__adaptive_random( std::vector< std::vector<int> > &neighb ) const
{
	int swarm_size = neighb.size();
	int pidx;		// for iterating over particles
	int nidx, j;	// for iterating over particles being connected to
	
	// clear previously defined topology
	for( pidx = 0; pidx < swarm_size; pidx++ )
		neighb[pidx].clear();
	
	// define new topology
	for( pidx = 0; pidx < swarm_size; pidx++ ){
		
		// the particle always connects back to itself, thus guaranteeing a minimum indegree of 1
		neighb[pidx].push_back( pidx );
		
		for( j = 1; j < m_neighb_param; j++ ){
			nidx = m_drng() * swarm_size;
			neighb[nidx].push_back( pidx );
			// No check performed to see whether pidx is already in neighb[nidx],
			// leading to a performance penalty in particle__get_best_neighbor() when it occurs.
		}
	}
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
	s << "topology:" << m_neighb_type << ' ';
	if( m_neighb_type == 2 || m_neighb_type == 4 )
		s << "topology param.:" << m_neighb_param << ' ';
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::pso)
