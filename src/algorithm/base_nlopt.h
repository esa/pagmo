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

#ifndef PAGMO_ALGORITHM_BASE_NLOPT_H
#define PAGMO_ALGORITHM_BASE_NLOPT_H

#include <cstddef>
#include <nlopt.h>
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
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE base_nlopt: public base
{
	protected:
		base_nlopt(nlopt_algorithm, bool, int, const double &);
		void evolve(population &) const;
		std::string human_readable_extra() const;
	private:
		struct nlopt_wrapper_data
		{
			problem::base const		*prob;
			decision_vector			*x;
			fitness_vector			*f;
			constraint_vector		*c;
			problem::base::c_size_type	c_comp;
		};
		static double objfun_wrapper(int, const double *, double *, void *);
		static double constraints_wrapper(int, const double *, double *, void *);
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & const_cast<nlopt_algorithm &>(m_algo);
			ar & const_cast<bool &>(m_constrained);
			ar & const_cast<std::size_t &>(m_max_iter);
			ar & const_cast<double &>(m_tol);
			ar & m_last_status;
		}  
		const nlopt_algorithm	m_algo;
		const bool		m_constrained;
		const std::size_t	m_max_iter;
		const double		m_tol;
		mutable int		m_last_status;
};

}}

#endif
