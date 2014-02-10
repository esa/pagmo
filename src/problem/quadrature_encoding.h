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

#ifndef PAGMO_PROBLEM_QUADRATURE_ENCODING_H
#define PAGMO_PROBLEM_QUADRATURE_ENCODING_H

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"

#define PI boost::math::constants::pi<double>()

namespace pagmo { namespace problem {

/**
 * Quadrature encoding problem. Transforms genes that encode angles to 
 * quadrature encoding. The resulting problem has two genes (i,j) for every
 * tranformed gene (x) of the original problem:
 * i = sin(x)
 * j = cos(x)
 * The sin component remains at the position of the original gene. The cos
 * component is added to the end of the chromosome.
 *
 * @author Daniel Hennes (daniel.hennes@gmail.com)
 */

class __PAGMO_VISIBLE quadrature_encoding : public base
{
public:

	/**
	 * Constructs the quadrature encoding on top of the original problem.
	 * @param[in] original problem
	 * @param[in] indices of genes in original chromosome to be transformed
	 */
	quadrature_encoding(const base &p = schwefel(1), const std::vector<size_type>&idx = std::vector<size_type>()):
		base(p.get_dimension() + idx.size(), p.get_i_dimension(), p.get_f_dimension(), p.get_c_dimension(), p.get_ic_dimension(), p.get_c_tol()), 
		m_original_problem(p.clone()), 
		m_idx(idx.begin(), idx.end())
	{
		assert(!m_idx.empty()); // we should have at least one index ...
		assert(*m_idx.rbegin() < m_original_problem->get_dimension()); // ... and largest index smaller than dimension of the problem

		// get bounds of original problem
		decision_vector lb = p.get_lb();
		decision_vector ub = p.get_ub();

		// set bounds of transformed problem for sin component
		std::set<size_type>::iterator it = m_idx.begin();
		for(size_type i = 0; i != lb.size(); i++)
		{
			if (*it == i)
			{
				lb[i] = -1.0;
				ub[i] = 1.0;
				it++;
			}
		}

		// .. for cos component
		for(size_type i = 0; i < m_idx.size(); i++)
		{
			lb.push_back(-1.0);
			ub.push_back(1.0);
		}

		// set new bounds
		set_bounds(lb, ub);
	}
	
	base_ptr clone() const { return base_ptr(new quadrature_encoding(*this)); }
	std::string get_name() const { return m_original_problem->get_name() + " [quadrature encoding]"; }

	/**
	 * Transforms decision vector from quadrature encoding to original encoding.
	 * @param[in] decision vector in quadrature encoding
	 * @param[out] decision vector in original encoding
	 */
	decision_vector transform2old(const decision_vector &x) const 
	{ 
		decision_vector _x(x);
		size_type ic = m_original_problem->get_dimension();
		for (std::set<size_type>::iterator it = m_idx.begin(); it != m_idx.end(); it++)
		{
			_x[*it] = atan2(x[*it], x[ic]);
			//_x[*it] = atan2(x[*it]*sqrt(1-.5*x[ic]*x[ic]), x[ic]*sqrt(1-.5*x[*it]*x[*it])); // this mapping produces a "more uniform" sampling
			ic++;
		}
		for (size_type i = 0; i < m_idx.size(); i++) _x.pop_back(); // remove unused genes
		return _x;
	}

	/**
	 * Transforms decision vector from original encoding to quadrature encoding.
	 * @param[in] decision vector in original encoding
	 * @param[out] decision vector in quadrature encoding
	 */
	decision_vector transform2new(const decision_vector &x) const 
	{ 
		decision_vector _x(x);
		for (std::set<size_type>::iterator it = m_idx.begin(); it != m_idx.end(); it++)
		{
			_x[*it] = sin(x[*it]);
			_x.push_back(cos(x[*it]));
		}
		return _x;
	}

protected:
	void objfun_impl(fitness_vector &f, const decision_vector &x) const 
	{
		m_original_problem->objfun(f, transform2old(x)); // will repeat safety checks, maybe use objfun_impl instead
	}

	base_ptr m_original_problem;
	std::set<size_type> m_idx;

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
	}

};
		
}} // namespaces

#endif // PAGMO_PROBLEM_QUADRATURE_ENCODING_H
