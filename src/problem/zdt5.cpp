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

#include <cmath>

#include "../exceptions.h"
#include "../types.h"
#include "../population.h"
#include "zdt5.h"

namespace pagmo { namespace problem {

/**
 * Will construct ZDT5.
 *
 * @param[in] bnr number of bitstrings used for the problem (integer dimension is 30 + 5 * (bnr-1))
 *
 * @see problem::base constructors.
 */
zdt5::zdt5(int bnr):base_unc_mo(30 + 5 * (bnr-1),30 + 5 * (bnr-1),2)
{
	// Set bounds.
	set_lb(0);
	set_ub(1);

		if (bnr <= 0) {
				pagmo_throw(value_error,"invalid dimension(s)");
		}

}

/// Clone method.
base_ptr zdt5::clone() const
{
	return base_ptr(new zdt5(*this));
}


/// Convergence metric for a decision_vector (0 = converged to the optimal front)
double zdt5::convergence_metric(const decision_vector &x) const
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


/// Implementation of the objective function.
void zdt5::objfun_impl(fitness_vector &f, const decision_vector &x) const
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

std::string zdt5::get_name() const
{
	return "ZDT5";
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::zdt5);
