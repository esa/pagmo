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

// 11/06/08 Created by Dario Izzo.

#include <cmath>
#include <sstream>
#include <vector>

#include "../../exceptions.h"
#include "base.h"
#include "classic.h"

//***********************************************************************************
//Classical problems
//***********************************************************************************

namespace pagmo
{
namespace problem {

test::test(int dim):base(std::vector<double>(dim, 0.0), std::vector<double>(dim, 1.0))
{
	if (dim <= 0) {
		pagmo_throw(value_error,"dimension of test problem must be positive");
	}
}

std::string test::id_object() const
{
	std::stringstream tmp;
	tmp << id_name() << "_" << getDimension();
	return tmp.str();
}

rastrigin::rastrigin(int dim):base(std::vector<double>(dim, -5.12), std::vector<double>(dim, 5.12))
{
	if (dim <= 0) {
		pagmo_throw(value_error,"dimension of rastrigin problem must be positive");
	}
}

std::string rastrigin::id_object() const
{
	std::stringstream tmp;
	tmp << id_name() << "_" << getDimension();
	return tmp.str();
}


schwefel::schwefel(int dim):base(std::vector<double>(dim, -500.0), std::vector<double>(dim, 500.0))
{
	if (dim <= 0) {
		pagmo_throw(value_error,"dimension of schwefel problem must be positive");
	}
}

std::string schwefel::id_object() const
{
	std::stringstream tmp;
	tmp << id_name() << "_" << getDimension();
	return tmp.str();
}


ackley::ackley(int dim):base(std::vector<double>(dim, -15.0), std::vector<double>(dim, 30.0))
{
	if (dim <= 0) {
		pagmo_throw(value_error,"dimension of ackley problem must be positive");
	}
}

std::string ackley::id_object() const
{
	std::stringstream tmp;
	tmp << id_name() << "_" << getDimension();
	return tmp.str();
}


rosenbrock::rosenbrock(int dim):base(std::vector<double>(dim, -5.0), std::vector<double>(dim, 10.0))
{
	if (dim <= 0) {
		pagmo_throw(value_error,"dimension of rosenbrock problem must be positive");
	}
}

std::string rosenbrock::id_object() const
{
	std::stringstream tmp;
	tmp << id_name() << "_" << getDimension();
	return tmp.str();
}


lennardjones::lennardjones(int atoms):base(3*atoms-6)
{
	if (atoms <= 0 || atoms < 3) {
		pagmo_throw(value_error,"number of atoms for lennard-jones problem must be greater than 2");
	}
	for (size_t i = 0; i < LB.size(); i++) {
		if ( (i != 0) && (i % 3) == 0 ) {
			LB[i] = 0.0;
			UB[i] = 6.0;
		} else {
			LB[i] = -3.0;
			UB[i] = 3.0;
		}
	}
}

std::string lennardjones::id_object() const
{
	std::stringstream tmp;
	tmp << id_name() << "_" << ((getDimension() + 6) / 3);
	return tmp.str();
}


levy::levy(int dim):base(std::vector<double>(dim, -10.0), std::vector<double>(dim, 10.0))
{
	if (dim <= 0) {
		pagmo_throw(value_error,"dimension of levy problem must be positive");
	}
}

std::string levy::id_object() const
{
	std::stringstream tmp;
	tmp << id_name() << "_" << getDimension();
	return tmp.str();
}

griewank::griewank(int dim):base(std::vector<double>(dim, -600.0), std::vector<double>(dim, 600.0))
{
	if (dim <= 0) {
		pagmo_throw(value_error,"dimension of griewank problem must be positive");
	}
}

std::string griewank::id_object() const
{
	std::stringstream tmp;
	tmp << id_name() << "_" << getDimension();
	return tmp.str();
}

}
}
