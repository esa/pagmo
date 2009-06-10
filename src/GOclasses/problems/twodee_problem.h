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

// 09/09/09 Created by Francesco Biscani.

#ifndef PAGMO_TWODEE_PROBLEM_H
#define PAGMO_TWODEE_PROBLEM_H

#include <string>
#include <vector>

#include "../../../config.h"
#include "GOproblem.h"

// Twodee problem.
class __PAGMO_VISIBLE twodee_problem: public GOProblem {
	public:
		twodee_problem(int);
		twodee_problem(int, const std::string &);
		virtual twodee_problem *clone() const {return new twodee_problem(*this);}
		virtual void post_evolution() const;
		virtual std::string id_object() const { return id_name(); }
	private:
		virtual double objfun_(const std::vector<double> &) const;
		mutable size_t		m_random_seed;
		const std::string	m_arguments;
};

#endif
