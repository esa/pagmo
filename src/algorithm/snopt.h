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

class __PAGMO_VISIBLE snopt: public base
{
public:
	snopt(const int major,const double feas=1e-10, const double opt = 1e-4);
	base_ptr clone() const;
	void evolve(population &) const;
protected:
	std::string human_readable_extra() const;

private:

	const int m_major;
	const double m_feas;
	const double m_opt;

	//This vector is here to allow the static function snopt_function_ to call the
	//GOProblem->objfun() without allocating memory for the chromosome. di_comodo will
	//be resized by the evolve command and then a refernece will be passed by snopt_a to
	//snopt_function_
	mutable decision_vector di_comodo;
};

}} //namespaces

#endif
