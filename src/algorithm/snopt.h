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

#ifndef PAGMO_ALGORITHM_SNOPT_H
#define PAGMO_ALGORITHM_SNOPT_H

#include "../config.h"
#include "../problem/base.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"

namespace pagmo { namespace algorithm {

/// Wrapper for the SNOPT solver
/**
 * SNOPT is a quite popular commercial solver coded in FORTRAN77 by GILL and MURRAY.
 * Its popularity stems from its being able to solve efficiently many different problems
 * thanks to what many would call 'black-magic' heuristic implemented in the solver that is
 * otherwise an SQP solver. We provide in PaGMO the wrappers around the libraries that
 * the user needs to have installed and licenced for in his computer.
 *
 * In order to interface to SNOPT succesfully, PaGMO needs to find in the system the libraries:
 * snopt, snprint, blas, f2c, m and gfortran (in case the snopt libraries were compiled using gfortran)
 *
 * CAREFUL:
 *	1 - SNOPT works only for minimization.
 *	2 - The final solution is guaranteed to be within the box constraints by forcing it after the snopt call
 *
 * From the SNOPT User-Manual:
 *
 * SNOPT is a general-purpose system for constrained optimization. It minimizes a
 * linear or nonlinear function subject to bounds on the variables and sparse linear or
 * nonlinear constraints. It is suitable for large-scale linear and quadratic programming
 * and for linearly constrained optimization, as well as for general nonlinear programs.
 * SNOPT finds solutions that are locally optimal, and ideally any nonlinear functions
 * should be smooth and users should provide gradients. It is often more widely useful.
 * For example, local optima are often global solutions, and discontinuities in the function
 * gradients can often be tolerated if they are not too close to an optimum. Unknown
 * gradients are estimated by finite differences.
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 *
 */

class __PAGMO_VISIBLE snopt: public base
{
public:

	snopt(const int major = 100,const double feas=1e-10, const double opt = 1e-4);
	base_ptr clone() const;
	void evolve(population &) const;
	void file_output(const bool);
	std::string get_name() const;

	//This structure contains one decision vector and one constraint vector as to allow
	//the static snopt function not to allocate any memory.
	struct preallocated_memory{
		decision_vector x;
		constraint_vector c;
		fitness_vector f;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & x;
			ar & c;
			ar & f;
		}
	};
protected:
	std::string human_readable_extra() const;

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<int &>(m_major);
		ar & const_cast<double &>(m_feas);
		ar & const_cast<double &>(m_opt);
		ar & m_file_out;
		ar & di_comodo;
	}  
	const int m_major;
	const double m_feas;
	const double m_opt;
	bool m_file_out;

	mutable preallocated_memory di_comodo;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::snopt)

#endif
