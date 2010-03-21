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

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../types.h"
#include <limits.h>
#include "base.h"
#include "ipopt.h"
#include "ipopt_cpp_wrapper/ipopt_problem.h"
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>


namespace pagmo { namespace algorithm {

ipopt::ipopt(const int &iter,const double &tol, const double &obj_tol) : m_iter(iter),m_tol(tol),m_obj_tol(obj_tol),m_screen_out(false)
{
	if (iter < 0) {
		pagmo_throw(value_error,"number of maximum iterations cannot be negative");
	}
	if (tol < 0 || tol > 1) {
		pagmo_throw(value_error,"tolerance is not in ]0,1[");
	}
	if (obj_tol < 0 || obj_tol > 1) {
		pagmo_throw(value_error,"obj_tol is not in ]0,1[");
	}
	app = IpoptApplicationFactory();

}

/// Clone method.
base_ptr ipopt::clone() const
{
	return base_ptr(new ipopt(*this));
}

/// Evolve implementation.
/**
 * Run IPOPT with the parameters specified in the constructor
 * At the end velocity is updated
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void ipopt::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type D = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension(), prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = pop.size();
	const problem::base::size_type Dc = D - prob_i_dimension;
	const std::vector<double>::size_type D_ineqc = prob.get_ic_dimension();
	const std::vector<double>::size_type D_eqc = prob_c_dimension - D_ineqc;
	const std::string name = prob.get_name();

	//We perform some checks to determine wether the problem/population are suitable for IPOPT
	if ( prob_i_dimension != 0  ) {
		pagmo_throw(value_error,"No integer part allowed yet....");
	}

	if ( Dc == 0  ) {
		pagmo_throw(value_error,"No continuous part....");
	}

	if ( prob_f_dimension != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and IPOPT is not suitable to solve it");
	}

	// Get out if there is nothing to do.
	if (NP == 0 || m_iter == 0) {
		return;
	}

	//create an instance of the ipopt_problem
	SmartPtr<Ipopt::TNLP> pagmo_nlp = new ipopt_problem(pop);

	// Change some options
	app->Options()->SetNumericValue("max_iter", m_iter);
	app->Options()->SetNumericValue("tol", m_tol);
	app->Options()->SetNumericValue("acceptable_obj_change_tol", m_obj_tol);

	// Intialize the IpoptApplication and process the options
	ApplicationReturnStatus status;
	status = app->Initialize();
	if (status != Solve_Succeeded) {
		pagmo_throw(value_error, "Error during IPOPT initialization!");
	}

	// Ask Ipopt to solve the problem
	status = app->OptimizeTNLP(pagmo_nlp);


	//Save the final point
	//for (integer i=0;i<n;i++) di_comodo.x[i] = x[i];
	//decision_vector newx = di_comodo.x;
	//std::transform(di_comodo.x.begin(), di_comodo.x.end(), pop.get_individual(bestidx).cur_x.begin(), di_comodo.x.begin(),std::minus<double>());
	//pop.set_x(bestidx,newx);
	//pop.set_v(bestidx,di_comodo.x);


}

/// Activate screen output
/**
 * Activate IPOPT screen output at a default level
 *
 * @param[in] p true or false
 */
void ipopt::screen_output(const bool p) {m_screen_out = p;}


std::string ipopt::human_readable_extra() const
{
	std::ostringstream s;
	s << "IPOPT - Max Iterations: " << m_iter << ", Tolerance: "<<m_tol<< ", Objective Function Tolerance: "<<m_obj_tol << std::endl;
	return s.str();
}

}} //namespaces
