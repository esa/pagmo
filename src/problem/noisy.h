/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
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

#ifndef PAGMO_PROBLEM_NOISY_H
#define PAGMO_PROBLEM_NOISY_H

#include <string>

#include "../serialization.h"
#include "ackley.h"
#include "../types.h"
#include "base.h"
#include "base_stochastic.h"

#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

namespace pagmo{ namespace problem {

/// Noisy meta-problem
/**
 * A meta-problem that transforms a problem into its stochastic
 * version by injecting noises to the fitness vector (and constraint vector)
 * with either uniform or gaussian distribution.
 *
 * @author Yung-Siang Liau (liauys@gmail.com)
 */
class __PAGMO_VISIBLE noisy : public base_stochastic
{
	public:
		/// Distribution type of the noise
		struct noise_distribution{
			enum type {NORMAL = 0, UNIFORM = 1};
		};
		//constructors
		noisy(const base & = ackley(1),
			  const double param_first = 0.0,
			  const double param_second = 0.1,
			  noise_distribution::type noise_type = noise_distribution::NORMAL,
			  unsigned int seed = 0);
		
		//copy constructor
		noisy(const noisy &);
		base_ptr clone() const;
		std::string get_name() const;

		void set_noise_param(double, double);
		double get_param_first();
		double get_param_second();

	protected:
		std::string human_readable_extra() const;
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void compute_constraints_impl(constraint_vector &, const decision_vector &) const;

	private:
		void inject_noise_f(fitness_vector&) const;
		void inject_noise_c(constraint_vector&) const;

		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_stochastic>(*this);
			ar & m_original_problem;
			ar & m_normal_dist;
			ar & m_uniform_dist;
			ar & m_param_first;
			ar & m_param_second;
			ar & m_noise_type;
		}

		base_ptr m_original_problem;
		mutable boost::normal_distribution<double> m_normal_dist;
		mutable boost::random::uniform_real_distribution<double> m_uniform_dist;
		double m_param_first;
		double m_param_second;
		noise_distribution::type m_noise_type;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::noisy);

#endif // PAGMO_PROBLEM_NOISY_H
