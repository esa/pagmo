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

namespace pagmo{ namespace problem {

/// Serializable boost::random distribution
/**
 * Implements a serializable boost::normal_distribution,
 * which could then be used safely within the PaGMO framework.
 *
 */
//TODO: If useful, can consider exposing this to a wider namespace, like rng_double?
class __PAGMO_VISIBLE normal_distribution_serializable: public boost::normal_distribution<double>
{
	friend class boost::serialization::access;
	public:
		/// Return value of the generator.
		typedef boost::normal_distribution<double>::input_type input_type;
		typedef boost::normal_distribution<double>::result_type result_type;
		/// Constructor from mean and sigma params
		/**
		 * Will invoke the corresponding base constructor.
		 */	
		normal_distribution_serializable(const double mean, const double sigma):boost::normal_distribution<double>(mean, sigma) {}

	private:
		// Serialization exploits the fact that the state of Boost distribution can be sent/received to/from standard streams.
		template <class Archive>
		void save(Archive &ar, const unsigned int) const
		{
			std::stringstream ss;
			ss << *static_cast<boost::normal_distribution<double> const *>(this);
			std::string tmp(ss.str());
			ar << tmp;
		}
		template <class Archive>
		void load(Archive &ar, const unsigned int)
		{
			std::string tmp;
			ar >> tmp;
			std::stringstream ss(tmp);
			ss >> *static_cast<boost::normal_distribution<double> *>(this);
		}
		BOOST_SERIALIZATION_SPLIT_MEMBER();
};


/// Noisy meta-problem
/**
 * Implements a meta-problem that transforms a problem into its stochastic 
 * version by injecting random noises to the fitness vector (and constraint vector)
 * with a certain distribution.
 *
 * @author Yung-Siang Liau (liauys@gmail.com)
 */
class __PAGMO_VISIBLE noisy : public base_stochastic
{
	public:
		//constructors
		noisy(const base & = ackley(1));
		noisy(const base &, const double mu, const double sigma, unsigned int seed);
		
		//copy constructor
		noisy(const noisy &);
		base_ptr clone() const;
		std::string get_name() const;

		void set_noise_param(double, double);
		double get_param_mu();
		double get_param_sigma();

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
			ar & m_param_mu;
			ar & m_param_sigma;
		}
		base_ptr m_original_problem;			
		mutable normal_distribution_serializable m_normal_dist;
		double m_param_mu;
		double m_param_sigma;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::noisy);

#endif // PAGMO_PROBLEM_NOISY_H
