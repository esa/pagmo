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
 * of the velocity update rule differing on the definition of the random vectors \f$r_1\f$ and \f$r_2\f$
 * \li Variant 1: \f$\mathbf r_1 = [r_1, r_1, ..., r_1]\f$, \f$\mathbf r_2 = [r_2, r_2, ..., r_2]\f$
 * \li Variant 2: \f$\mathbf r_1 = [r_1, r_1, ..., r_1]\f$, \f$\mathbf r_2 = [r_1, r_1, ..., r_1]\f$
 * \li Variant 3: \f$\mathbf r_1 = [r_{11}, r_{12}, ..., r_{1n}]\f$, \f$\mathbf r_2 = [r_{21}, r_{21}, ..., r_{2n}]\f$
 * \li Variant 4: \f$\mathbf r_1 = [r_{11}, r_{12}, ..., r_{1n}]\f$, \f$\mathbf r_2 = [r_{11}, r_{11}, ..., r_{1n}]\f$
 *
 * At each call of the evolve method a number of function evaluations equal to m_gen * pop.size()
 * is performed.
 *
 * The algorithm is suitable for box-constrained single-objective continuous optimization.
 *
 * @see http://swarmintelligence.org/ for a fair descritpion of the algorithm
 * @see http://www.engr.iupui.edu/~shi/Coference/psopap4.html for the first paper on this algorithm
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 */

class __PAGMO_VISIBLE pso: public base
{
public:
	pso(int gen = 1, double omega = 0.65, double eta1 = 2.0, double eta2 = 2.0, double vcoeff = 0.2, int variant = 3);
	base_ptr clone() const;
	void evolve(population &) const;
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
	}  
	// Number of generations.
	const int m_gen;
	// Particle Inertia
	const double m_omega;
	// Weight of the social component
	const double m_eta1;
	// Weight of the social component
	const double m_eta2;
	// Velocity coefficient
	const double m_vcoeff;
	// Startegy
	const int m_variant;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::pso);

#endif // PSO_H
