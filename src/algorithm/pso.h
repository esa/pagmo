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

#ifndef PAGMO_ALGORITHM_PSO_H
#define PAGMO_ALGORITHM_PSO_H

#include "../config.h"
#include "../serialization.h"
#include "base.h"


namespace pagmo { namespace algorithm {

/// Particle Swarm optimization
/**
 *
 * Particle swarm optimization (PSO) is a population based algorithm that has been proposed in the
 * mid nineties and that is inspired by the foraging behaviour of swarms. In PSO each point has
 * memory of the position where it achieved the best performance \f$\mathbf x^l_i\f$ and of the swarm
 * 'champion' position \f$ \mathbf x^g \f$ and uses this information to update its position using the equation:
 * \f[
 *	\mathbf v_{i+1} = \omega \mathbf v_i + \eta_1 \mathbf r_1 \cdot \left( \mathbf x_i - \mathbf x^l_i \right)
 *	+ \eta_2 \mathbf r_2 \cdot \left(  \mathbf x_i - \mathbf x^g \right)
 * \f]
 * \f[
 *	\mathbf x_{i+1} = \mathbf x_i + \mathbf v_i
 * \f]
 *
 * The user can specify the values for \f$\omega, \eta_1, \eta_2\f$ and the magnitude of the maximum velocity
 * allowed. this last value is evaluated for each search direction as the product of \f$ vcoeff\f$ and the
 * search space width along that direction. The user can also specify one of four variants:
 *
 * the velocity update rule differing on the definition of the random vectors \f$r_1\f$ and \f$r_2\f$
 * \li Variant 1: \f$\mathbf r_1 = [r_{11}, r_{12}, ..., r_{1n}]\f$, \f$\mathbf r_2 = [r_{21}, r_{21}, ..., r_{2n}]\f$
 * \li Variant 2: \f$\mathbf r_1 = [r_{11}, r_{12}, ..., r_{1n}]\f$, \f$\mathbf r_2 = [r_{11}, r_{11}, ..., r_{1n}]\f$
 * \li Variant 3: \f$\mathbf r_1 = [r_1, r_1, ..., r_1]\f$, \f$\mathbf r_2 = [r_2, r_2, ..., r_2]\f$
 * \li Variant 4: \f$\mathbf r_1 = [r_1, r_1, ..., r_1]\f$, \f$\mathbf r_2 = [r_1, r_1, ..., r_1]\f$
 *
 * \li Variant 5: \f$\mathbf r_1 = [r_1, r_1, ..., r_1]\f$, \f$\mathbf r_2 = [r_1, r_1, ..., r_1]\f$
 *
 * At each call of the evolve method a number of function evaluations equal to m_gen * pop.size()
 * is performed.
 * 
 * The algorithm is suitable for box-constrained single-objective continuous optimization.
 * 
 * @see http://www.particleswarm.info/ for a repository of information related to PSO
 * @see http://dx.doi.org/10.1007/s11721-007-0002-0 for a recent survey
 * @see http://www.engr.iupui.edu/~shi/Coference/psopap4.html for the first paper on this algorithm
 * 
 * @author Dario Izzo (dario.izzo@googlemail.com)
 * @author Luis Simoes (luis.f.m.simoes@gmail.com)
 */

class __PAGMO_VISIBLE pso: public base
{
public:
	pso(int gen=1, double omega = 0.7298, double eta1 = 2.05, double eta2 = 2.05, double vcoeff = 0.5, int variant = 5, int neighb_type = 2, int neighb_param = 4 );
	base_ptr clone() const;
	void evolve(population &) const;
	decision_vector particle__get_best_neighbor( population::size_type pidx, std::vector< std::vector<int> > &neighb, const std::vector<decision_vector> &lbX, const std::vector<fitness_vector> &lbfit, const problem::base &prob ) const;
	void initialize_topology__gbest( const population &pop, decision_vector &gbX, fitness_vector &gbfit, std::vector< std::vector<int> > &neighb ) const;
	void initialize_topology__lbest( std::vector< std::vector<int> > &neighb ) const;
	void initialize_topology__von( std::vector< std::vector<int> > &neighb ) const;
	void initialize_topology__adaptive_random( std::vector< std::vector<int> > &neighb ) const;
	std::string get_name() const;
protected:
	std::string human_readable_extra() const;
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<int &>(m_gen);
		ar & const_cast<double &>(m_omega);
		ar & const_cast<double &>(m_eta1);
		ar & const_cast<double &>(m_eta2);
		ar & const_cast<double &>(m_vcoeff);
		ar & const_cast<int &>(m_variant);
		ar & const_cast<int &>(m_neighb_type);
		ar & const_cast<int &>(m_neighb_param);
	}  
	// Number of generations
	const int m_gen;
	// Particle Inertia weight, or alternatively the constriction coefficient
	const double m_omega;
	// magnitude of the force, applied to the particle's velocity, in the direction of its previous best position
	const double m_eta1;
	// magnitude of the force, applied to the particle's velocity, in the direction of the best position in its neighborhood
	const double m_eta2;
	// Velocity coefficient: velocity values will range in [ - var_range * m_vcoeff, var_range * m_vcoeff ]
	const double m_vcoeff;
	// Velocity update formula
	const int m_variant;
	// Swarm topology
	const int m_neighb_type;
	// parameterization of the swarm topology
	const int m_neighb_param;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::pso)

#endif // PSO_H
