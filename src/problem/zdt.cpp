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

#include <cmath>
#include <boost/math/constants/constants.hpp>

#include "../exceptions.h"
#include "../types.h"
#include "../population.h"
#include "zdt.h"

static int __check__(int N)
{
		if (N > 6 || N < 1) {
				pagmo_throw(value_error, "the problem id needs to be one of [1..6]");
		}
		return N;
}

namespace pagmo { namespace problem {

/**
 * Will construct a problem of the ZDT test-suite.
 * @param[in] id problem number: 1-6 possible.
 * @param[in] param_1 problem parameter, dimension for ZDT1,2,3,4,6 - number of binary strings for ZDT5
 *
 * @see problem::base constructors.
 */
zdt::zdt(size_type id, size_type param_1)
	:base_unc_mo( (__check__(id) != 5) ? param_1 : 30 + 5 * (param_1 - 1), (__check__(id) != 5) ? 0 : 30 + 5 * (param_1 - 1), 2 ),
	m_problem_number(__check__(id))
{
	if(param_1 <= 0) {
		pagmo_throw(value_error,"Invalid parameter (problem dimension needs to be higher)");
	}

		// set the bounds for the current problem
		switch(m_problem_number)
		{
		case 1:
		case 2:
		case 3:
		case 6:
		{
			set_lb(0.0);
			set_ub(1.0);
			break;
		}
		case 4:
		{
			set_lb(-5.0);
			set_ub(5.0);
			set_lb(0,0.0);
			set_ub(0,1.0);
			break;
		}
		case 5:
		{
			set_lb(0);
			set_ub(1);
			break;
		}
		default:
			pagmo_throw(value_error, "Error: There are only 6 test functions in this test suite!");
			break;
		}
}

/// Clone method.
base_ptr zdt::clone() const
{
	return base_ptr(new zdt(*this));
}

/// Convergence metric for a decision_vector (0 = converged to the optimal front)
double zdt::convergence_metric(const decision_vector &x) const
{
	switch(m_problem_number)
	{
	case 1:
	case 2:
	case 3:
		return g0123_convergence_metric(x);
	case 4:
		return g04_convergence_metric(x);
	case 5:
		return g05_convergence_metric(x);
	case 6:
		return g06_convergence_metric(x);
	default:
		pagmo_throw(value_error, "Error: There are only 6 test functions in this test suite!");
	}
	return -1.0;
}

/// Implementation of the objective function.
void zdt::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	switch(m_problem_number)
	{
	case 1:
		g01_objfun_impl(f,x);
		break;
	case 2:
		g02_objfun_impl(f,x);
		break;
	case 3:
		g03_objfun_impl(f,x);
		break;
	case 4:
		g04_objfun_impl(f,x);
		break;
	case 5:
		g05_objfun_impl(f,x);
		break;
	case 6:
		g06_objfun_impl(f,x);
		break;
	default:
		pagmo_throw(value_error, "Error: There are only 6 test functions in this test suite!");
		break;
	}
}

// shared convergence metric for ZDT1, ZDT2 and ZDT3
double zdt::g0123_convergence_metric(const decision_vector &x) const
{
	double c = 0.0;
	double g = 0.0;

	for(problem::base::size_type j = 1; j < x.size(); ++j) {
		g += x[j];
	}
	c += 1 + (9 * g) / (x.size()-1);
	return c - 1;
}

void zdt::g01_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 2);
	pagmo_assert(x.size() == get_dimension());

	double g = 0;

	f[0] = x[0];

	for(problem::base::size_type i = 1; i < x.size(); ++i) {
		g += x[i];
	}
	g = 1 + (9 * g) / (x.size()-1);

	f[1] = g * ( 1 - sqrt(x[0]/g));

}

void zdt::g02_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
		pagmo_assert(f.size() == 2);
		pagmo_assert(x.size() == get_dimension());

		double g = 0;

		f[0] = x[0];

		for(problem::base::size_type i = 1; i < x.size(); ++i) {
				g += x[i];
		}
		g = 1 + (9 * g) / (x.size() - 1);

		f[1] = g * ( 1 - (x[0]/g)*(x[0]/g));
}

void zdt::g03_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
		pagmo_assert(f.size() == 2);
		pagmo_assert(x.size() == get_dimension());

		double g = 0;

		f[0] = x[0];

		for(problem::base::size_type i = 1; i < x.size(); ++i) {
				g += x[i];
		}
		g = 1 + (9 * g) / (x.size()-1);

		f[1] = g * ( 1 - sqrt(x[0]/g) - x[0]/g * sin(10 * boost::math::constants::pi<double>() * x[0]));

}

double zdt::g04_convergence_metric(const decision_vector &x) const
{
		double c = 0.0;
		double g = 0.0;

		for(problem::base::size_type j = 1; j < x.size(); ++j) {
				g += x[j]*x[j] - 10 * cos(4 * boost::math::constants::pi<double>() * x[j]);
		}
		c += 1 + 10 * (x.size()-1) + g;
		return c  - 1;
}

// ZDT4': lambda x: 1 + 10 * (len(x) - 1) + sum([xi**2 - 10*np.cos(4*np.pi*xi) for xi in x[1:]])
void zdt::g04_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
		pagmo_assert(f.size() == 2);
	pagmo_assert(x.size() == get_dimension());

		double g = 1 + 10 * (x.size() - 1);

		f[0] = x[0];

		for(problem::base::size_type i = 1; i < x.size(); ++i) {
				g += x[i]*x[i] - 10 * cos(4 * boost::math::constants::pi<double>() * x[i]);
		}

		f[1] = g * ( 1 - sqrt(x[0]/g) );

}


double zdt::g05_convergence_metric(const decision_vector &x) const
{
		double c = 0.0;
		double g = 0.0;
		int k = 30;

		problem::base::size_type n_vectors = ((x.size()-30)/5)  +  1;
		std::vector<int> u(n_vectors,0);
		std::vector<int> v(n_vectors);

		for (problem::base::size_type i=1; i<n_vectors; i++)
		{
				for (int j=0; j<5; j++)
				{
						if (x[k] == 1)
						{
								 u[i]++;
						}
						k++;
				}
		}

		for (problem::base::size_type i=1; i<n_vectors; i++)
		{
				if (u[i] < 5)
				{
						v[i] = 2 + u[i];
				}
				else
				{
						v[i] = 1;
				}
		}

		for (problem::base::size_type i=1; i<n_vectors; i++)
		{
				g += v[i];
		}
		c +=  g;
		return c - n_vectors + 1;
}

void zdt::g05_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
		pagmo_assert(f.size() == 2);
				pagmo_assert(x.size() == get_dimension());
		double g = 0;
				problem::base::size_type size_x = x.size();
				problem::base::size_type n_vectors = ((size_x-30)/5)  +  1;
				int j, k = 30;
				std::vector<int> u(n_vectors);
				std::vector<int> v(n_vectors);


				for (problem::base::size_type i=0; i<n_vectors; i++)
				{
				   u[i] = 0;
				}
				for (j=0; j<30; j++)
				{
				  if (x[j] == 1)
				  {
						 u[0]++;
				  }
				}
				for (problem::base::size_type i=1; i<n_vectors; i++)
				{
				  for (j=0; j<5; j++)
				  {
						if (x[k] == 1)
						{
								u[i]++;
						}
						k++;
				  }
				}
				f[0] = 1.0 + u[0];
				for (problem::base::size_type i=1; i<n_vectors; i++)
				{
				   if (u[i] < 5)
				   {
						  v[i] = 2 + u[i];
				   }
				   else
				   {
						  v[i] = 1;
				   }
				}
				g = 0;
				for (problem::base::size_type i=1; i<n_vectors; i++)
				{
				   g += v[i];
				}

				f[1] = g * (1.0/f[0]);


}

double zdt::g06_convergence_metric(const decision_vector &x) const
{
		double c = 0.0;
		double g = 0.0;

		for(problem::base::size_type j = 1; j < x.size(); ++j) {
				g += x[j];
		}
		c += 1 + 9 * pow((g / (x.size()-1)),0.25);

		return c - 1;
}


void zdt::g06_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
		pagmo_assert(f.size() == 2);
		pagmo_assert(x.size() == get_dimension());

		double g = 0;

		f[0] = 1 - exp(-4*x[0])*pow(sin(6*boost::math::constants::pi<double>()*x[0]),6);

		for(problem::base::size_type i = 1; i < x.size(); ++i) {
				g += x[i];
		}
		g = 1 + 9 * pow((g/(x.size()-1)),0.25);

		f[1] = g * ( 1 - (f[0]/g)*(f[0]/g));

}


std::string zdt::get_name() const
{
	std::string retval("ZDT");
	retval.append(boost::lexical_cast<std::string>(m_problem_number));

	return retval;
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::zdt)
