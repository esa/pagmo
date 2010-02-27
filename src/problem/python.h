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

#ifndef PAGMO_PROBLEM_PYTHON_H
#define PAGMO_PROBLEM_PYTHON_H

#include <string>

#include "../config.h"
#include "../types.h"
#include "base.h"

namespace pagmo { namespace problem {

class __PAGMO_VISIBLE python: public base
{
	public:
		python(int, int, int, int, int);
		python(const double &, const double &, int, int, int, int, int);
		python(const decision_vector &, const decision_vector &, int, int, int, int);
		virtual fitness_vector py_objfun(const decision_vector &) const = 0;
		virtual std::string human_readable_extra() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
};

} }

#endif
