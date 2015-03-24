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

#ifndef PAGMO_PROBLEM_ROBUST_H
#define PAGMO_PROBLEM_ROBUST_H

#include <string>
#include <boost/functional/hash.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "../serialization.h"
#include "ackley.h"
#include "../types.h"
#include "base_stochastic.h"

namespace pagmo{ namespace problem {

/// Robust meta-problem
/**
 * Transforms any problem into a stochastic problem where the objective
 * function is evaluated sampling at random in a neighbourhood of the input
 * chromosome. The solution to the resulting problem is robust
 * to input noises in the given neighbourhood.
 *
 * @author Yung-Siang Liau (liauys@gmail.com)
 * @author Dario Izzo (dario.izzo@gmail.com)
 *
 */
class __PAGMO_VISIBLE robust : public base_stochastic
{
	public:

		//constructors
		robust(const base & = ackley(1),
			   unsigned int trials = 1,
			   const double param_rho = 0.1,
			   unsigned int seed = 0u);
		
		//copy constructor
		robust(const robust &);
		base_ptr clone() const;
		std::string get_name() const;

		void set_rho(double);
		double get_rho() const;

	protected:
		std::string human_readable_extra() const;
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void compute_constraints_impl(constraint_vector &, const decision_vector &) const;

	private:
		void inject_noise_x(decision_vector &) const;

		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_stochastic>(*this);
			ar & m_original_problem;
			ar & m_normal_dist;
			ar & m_uniform_dist;
			ar & m_trials;
			ar & m_rho;
		}

		base_ptr m_original_problem;
		mutable boost::normal_distribution<double> m_normal_dist;
		mutable boost::random::uniform_real_distribution<double> m_uniform_dist;
		unsigned int m_trials;
		double m_rho;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::robust)

#endif // PAGMO_PROBLEM_ROBUST_H
