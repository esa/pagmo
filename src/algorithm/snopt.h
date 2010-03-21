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

#ifndef PAGMO_ALGORITHM_SNOPT_H
#define PAGMO_ALGORITHM_SNOPT_H


#include "../config.h"
#include "base.h"
#include "../problem/base.h"
#include "snopt_cpp_wrapper/snopt_PAGMO.h"
#include "snopt_cpp_wrapper/snfilewrapper_PAGMO.h"

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
 * CAREFUL: SNOPT works only for minimization
 *
 * From the SNOPT User-Manual:
 *
 * SNOPT is a general-purpose system for constrained optimization. It minimizes a
 * linear or nonlinear function sub ject to bounds on the variables and sparse linear or
 * nonlinear constraints. It is suitable for large-scale linear and quadratic programming
 * and for linearly constrained optimization, as well as for general nonlinear programs.
 * SNOPT ﬁnds solutions that are locally optimal, and ideally any nonlinear functions
 * should be smooth and users should provide gradients. It is often more widely useful.
 * For example, local optima are often global solutions, and discontinuities in the function
 * gradients can often be tolerated if they are not too close to an optimum. Unknown
 * gradients are estimated by ﬁnite diﬀerences.
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 *
 */

class __PAGMO_VISIBLE snopt: public base
{
public:

	snopt(const int major,const double feas=1e-10, const double opt = 1e-4);
	base_ptr clone() const;
	void evolve(population &) const;
	void screen_output(const bool);
	void file_output(const bool);

	//This structure contains one decision vector and one constraint vector as to allow
	//the static snopt function not to allocate any memory.
	struct preallocated_memory{
		decision_vector x;
		constraint_vector c;
		fitness_vector f;
	};
protected:
	std::string human_readable_extra() const;

private:
	const int m_major;
	const double m_feas;
	const double m_opt;
	bool m_screen_out;
	bool m_file_out;

	mutable preallocated_memory di_comodo;
};

}} //namespaces

#endif
