/***************************************************************************
 *   Copyright (C) 2007, 2008 by Francesco Biscani   *
 *   bluescarni@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/pure_virtual.hpp>
#include <boost/utility.hpp>

#include "../src/GOclasses/basic/individual.h"
#include "../src/GOclasses/problems/TrajectoryProblems.h"

using namespace boost::python;

struct GOProblemWrap: GOProblem, wrapper<GOProblem>
{
	double objfun(const std::vector<double> &x)
	{
		return this->get_override("objfun")(x);
	}
};

// Instantiate the PyGMO module.
BOOST_PYTHON_MODULE(_PyGMO)
{
	class_<GOProblemWrap, boost::noncopyable> class_gop("goproblem", "Base GO problem", no_init);
	class_gop.def("objfun", pure_virtual(&GOProblem::objfun));
	class_gop.add_property("dimension", &GOProblem::getDimension, "Dimension of the problem.");

	class_<messengerfullProb, bases<GOProblem> > class_mfp("messenger_full_problem", "Messenger full problem.", init<>());

	class_<Individual> class_ind("individual", "Individual.", init<GOProblem &>());
	class_ind.add_property("fitness", &Individual::getFitness, "Fitness.");
}
