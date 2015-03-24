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

#ifndef PAGMO_SPHERES_H
#define PAGMO_SPHERES_H

#include <string>
#include <vector>
#include <gsl/gsl_odeiv2.h>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base_stochastic.h"
#include "../rng.h"

namespace pagmo { namespace problem {

/// Evolutionary Neuro-Controller for the MIT Spheres (perception-action defined in the absolute frame)
/**
 *
 * \image html spheres.jpg "The MIT spheres test-bed on board of the ISS."
 * \image latex spheres.jpg "The MIT spheres test-bed on boeard of the ISS." width=5cm
 *
 * This problem (a stochastic optimization problem) aims at 'evolving' an artificial neural network
 * able to steer the positions of three satellites representing the SPHERES MIT test-bed
 * on-board of the ISS. It requies GSL libraries and thus it is compiled only if that
 * option is activated when compiling pagmo.
 *
 * The objective function is the average over multiple runs of the following fitness
 * \f[
 * 	F = |L1^2 - r_1^2| + |L2^2 - r_2^2| + |L3^2 - r_3^2| + (|v_{f1}^2| + |v_{f2}^2| + |v_{f3}^2|) / 2
 * \f]
 *
 * rewarding neurocontrollers that drive the spheres towards a triangular configuration in space
 * with final zero absolute velocity. L1<L2<L3 are the three sides of the trianlgular formation desired and
 * \f$ r_1<r_2<r_3 \f$ the three actual interspheres distances at the end of one simulation.
 *
 * NOTE: the dynamical model of the spheres is here that of point masses. As a consequence, each sphere
 * perception and action (sensing and actuating) must be defined with respect to a reference frame that is known
 * by all spheres (a star-tracker would give this information). This seemingly small detail brings
 * to an important bias in learning (i.e. the final triangle is always achieved with the same absolute
 * orientation!!!). In pagmo::problem::spheres_q such a bias is removed by defining perception and action
 * in the sphere's body frame.
 *
 * @author Dario Izzo (dario.izzo@esa.int)
 */

class __PAGMO_VISIBLE spheres: public base_stochastic
{
	static int ode_func( double t, const double y[], double f[], void *params );
	public:
		/// Constructor
		/**
		 * Initializes the Sphere simulator, the neural net and everything else!!!
		 *
		 * @param[in] n_evaluations Number of initial conditions each neural network fitness is evaluated upon
		 * @param[in] n_hidden number of hidden neurons in the neural net
		 * @param[in] ode_prec precision requested to adapt the ode-solver step size
		 * @param[in] seed seed used to produce all random initial conditions
		 * @param[in] symmetric a boolean value that, if true, indicates that the neural network
		 * does not distinguish among permutations of its input values due to sphere ID exchange.
		 * @param[in] sim_time Time after wich the fitness is evaluated in the simualtion
		 * @param[in] sides The three sides of the trianglular formation to acquire and maintain

*/
		spheres(int n_evaluations = 10, int n_hidden = 10, double ode_prec = 1E-6, unsigned int seed = 0, bool symmetric = false, double sim_time = 50.0, const std::vector<double>& sides = std::vector<double>(3,0.5));

		/// Copy Constructor
		/**
		 * Necessary to properly handle the gsl variables (pointers!!!)
		 */
		spheres(const spheres &);

		/// Destructor
		/**
		 * Necessary to properly clear memory allocated by gsl routines
		 */
		~spheres();

		/// Post evaluation of the neural controller
		/**
		 * It tests a given neural controller over a large set of initial conditions (different from those
		 * the controller's fitness was actually evaluated upon during evolution)
		 *
		 * @param[in] x chromosome encoding the neural network
		 * @param[in] N number of initial conditions we want it to be tested against
		 * @param[in] seed seed used to produce all random initial conditions
		 * @return A vector containing, for each simulation, a vector of [ic, fitness] sorted by fitness
		 */
		std::vector<std::vector<double> > post_evaluate(const decision_vector &x, int N = 25000, unsigned int seed = 0) const;

		/// Performs one "detailed" simulation
		/**
		 * One single run of the ode solver is made, maintaining the states along the simulation ....)
		 *
		 * @param[in] x chromosome encoding the neural network
		 * @param[in] ic initial condition for the simulation
		 * @param[in] N number of points to retain along the simulation (plus the ic, of course)
		 * @return A vector containing, for each time, a vector of [t, states]
		 */
		std::vector<std::vector<double> > simulate(const decision_vector & x, const std::vector<double> &ic, int N) const;
		std::string get_name() const;
		base_ptr clone() const;
		
		/// Gets the weights of the neural network
		std::vector<double> get_nn_weights(decision_vector x) const;

	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		std::string human_readable_extra() const;
	private:
		// Class representing a feed forward neural network
		class ffnn {
				friend class spheres;
			public:
				ffnn(const unsigned int, const unsigned int,const unsigned int);
				void eval(double[], const double[]) const;
				void set_weights(const std::vector<double> &);
			private:
				friend class boost::serialization::access;
				template <class Archive>
				void serialize(Archive &ar, const unsigned int)
				{
					ar & const_cast<unsigned int &>(m_n_inputs);
					ar & const_cast<unsigned int &>(m_n_hidden);
					ar & const_cast<unsigned int &>(m_n_outputs);
					ar & m_weights;
					ar & m_hidden;
				}
				const unsigned int m_n_inputs;
				const unsigned int m_n_hidden;
				const unsigned int m_n_outputs;
				std::vector<double> m_weights;
				mutable std::vector<double> m_hidden;
		};
		void set_nn_weights(const decision_vector& x) const;
		double single_fitness( const std::vector<double> &, const ffnn& ) const;
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_stochastic>(*this);
			ar & m_ffnn;
			ar & m_n_evaluations;
			ar & m_n_hidden_neurons;
			ar & const_cast<double &>(m_numerical_precision);
			ar & m_ic;
			ar & m_symm;
			ar & m_sim_time;
			ar & m_sides;
		}
		gsl_odeiv2_driver*				m_gsl_drv_pntr;
		gsl_odeiv2_system				m_sys;
		mutable ffnn					m_ffnn;
		int								m_n_evaluations;
		int								m_n_hidden_neurons;
		const double					m_numerical_precision;
		mutable std::vector<double>		m_ic;	
		bool							m_symm;
		double							m_sim_time;
		std::vector<double>				m_sides;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::spheres)
#endif // PAGMO_SPHERES_H
