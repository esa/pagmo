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

#include <string>

#include "../types.h"
#include "base.h"
#include "python.h"

namespace pagmo { namespace problem {

python::python(int n, int ni, int nf, int nc, int nic):base(n,ni,nf,nc,nic) {}

python::python(const double &lb, const double &ub, int n, int ni, int nf, int nc, int nic):base(lb,ub,n,ni,nf,nc,nic) {}

python::python(const decision_vector &lb, const decision_vector &ub, int ni, int nf, int nc, int nic):base(lb,ub,ni,nf,nc,nic) {}

void python::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	f = py_objfun(x);
}

std::string python::human_readable_extra() const
{
	return base::human_readable_extra();
}

} }
