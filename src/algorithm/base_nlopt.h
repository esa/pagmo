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

#ifndef PAGMO_ALGORITHM_BASE_NLOPT_H
#define PAGMO_ALGORITHM_BASE_NLOPT_H

#include <cstddef>
#include <nlopt.hpp>
#include <string>

#include "../config.h"
#include "../population.h"
#include "../problem/base.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"

namespace pagmo { namespace algorithm {

/// Base class for wrapping NLopt's algorithms.
/**
 * NLopt is a free/open-source library for nonlinear optimization, providing a common interface for a number of different free optimization
 * routines available online as well as original implementations of various other algorithms.
 *
 * This class provides a common interface for NLopt algorithms. The type of algorithm that this class will wrap and the termination conditions
 * can be specified in the constructor.
 *
 * The evolve() method will select the best individual from the population and will optimise its continuous part using the specified NLopt algorithm.
 * After the optimisation, the individual will be unconditionally re-inserted in the same position in the population.
 *
 * All algorithms provided in NLopt are single-objective continuous minimisers.
 *
 * @see http://ab-initio.mit.edu/wiki/index.php/NLopt
 *
 * @author Francesco Biscani (bluescarni@gmail.com), Dario Izzo(dario.izzo@googlemail.com)
 */
class __PAGMO_VISIBLE base_nlopt: public base
{
	protected:
		base_nlopt(nlopt::algorithm, bool, bool, int, const double &, const double &);
		void evolve(population &) const;
		std::string human_readable_extra() const;
	private:
		struct nlopt_wrapper_data
		{
			problem::base const		*prob;
			decision_vector			x;
			decision_vector			dx;
			fitness_vector			f;
			constraint_vector		c;
			problem::base::c_size_type	c_comp;
		};
		int get_last_status() const;
		static double objfun_wrapper(const std::vector<double> &, std::vector<double> &, void*);
		static double constraints_wrapper(const std::vector<double> &, std::vector<double> &, void*);
		virtual void set_local(size_t d) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & const_cast<nlopt::algorithm &>(m_algo);
			ar & const_cast<bool &>(m_constrained);
			ar & const_cast<bool &>(m_only_ineq);
			ar & const_cast<std::size_t &>(m_max_iter);
			ar & const_cast<double &>(m_ftol);
			ar & const_cast<double &>(m_xtol);
		}
		const nlopt::algorithm	m_algo;
	protected:
		/// NLOPT optimization method
		mutable nlopt::opt	m_opt;
	private:
		const bool		m_constrained;
		const bool		m_only_ineq;
	protected:
		/// Maximum number of iterations
		const std::size_t	m_max_iter;
		/// Tolerance on the fitness function variation (stopping criteria)
		const double		m_ftol;
		/// Tolerance on the decision_vector variation function (stopping criteria)
		const double		m_xtol;
};

}}

#endif
