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

#include "ipopt.h"
#include "../exceptions.h"
#include "../types.h"


#include "ipopt_cpp_wrapper/ipopt_problem.h"
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>
#include <limits.h>


namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify some of the parameters of the IPOPT solver.
 *
 * @param[in] max_iter Maximum number of major iterations (refer to IPOPT manual).
 * @param[in] tol Desired convergence tolerance (relative). (refer to IPOPT manual).
 * @param[in] opt "Acceptance" stopping criterion based on objective function change. If the relative
 * change of the objective function (scaled by Max(1,|f(x)|)) is less than this value,
 * this part of the acceptable tolerance termination is satisfied. (refer to IPOPT manual).
 * @throws value_error if max_iter is not positive, and tol,acceptable_obj_change_tol are not in \f$]0,1[\f$
 */


ipopt::ipopt(const int &max_iter,const double &tol, const double &acceptable_obj_change_tol) : m_max_iter(max_iter),m_tol(tol),m_acceptable_obj_change_tol(acceptable_obj_change_tol),m_screen_out(false)
{
	if (max_iter < 0) {
		pagmo_throw(value_error,"number of maximum iterations cannot be negative");
	}
	if (tol < 0 || tol > 1) {
		pagmo_throw(value_error,"tolerance is not in ]0,1[");
	}
	if (acceptable_obj_change_tol < 0 || acceptable_obj_change_tol > 1) {
		pagmo_throw(value_error,"obj_tol is not in ]0,1[");
	}

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
	const problem::base::size_type D = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension(), prob_f_dimension = prob.get_f_dimension();
	const population::size_type NP = pop.size();
	const problem::base::size_type Dc = D - prob_i_dimension;
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
	if (NP == 0 || m_max_iter == 0) {
		return;
	}

	//create an instance of the ipopt_problem
	::Ipopt::SmartPtr< ::Ipopt::TNLP> pagmo_nlp = new ipopt_problem(&pop);

	//create an instance of the IpoptApplication
	::Ipopt::SmartPtr< ::Ipopt::IpoptApplication> m_app = new ::Ipopt::IpoptApplication();


	m_app->Options()->SetStringValue("hessian_approximation", "limited-memory");
	m_app->Options()->SetIntegerValue("print_level", 5);

	// Termination Criteria
	m_app->Options()->SetIntegerValue("max_iter", m_max_iter);
	m_app->Options()->SetNumericValue("tol", m_tol);
	m_app->Options()->SetNumericValue("acceptable_obj_change_tol", m_acceptable_obj_change_tol);

	// Intialize the IpoptApplication and process the options
	Ipopt::ApplicationReturnStatus status;
	status = m_app->Initialize();
	if (status != Ipopt::Solve_Succeeded) {
		pagmo_throw(value_error, "Error during IPOPT initialization!");
	}

	// Ask Ipopt to solve the problem
	status = m_app->OptimizeTNLP(pagmo_nlp);


	//Save the final point
	//for (integer i=0;i<n;i++) di_comodo.x[i] = x[i];
	//decision_vector newx = di_comodo.x;
	//std::transform(di_comodo.x.begin(), di_comodo.x.end(), pop.get_individual(bestidx).cur_x.begin(), di_comodo.x.begin(),std::minus<double>());
	//pop.set_x(bestidx,newx);
	//pop.set_v(bestidx,di_comodo.x);


}

/// Activate screen output
/**
 * Activate/Deactivate IPOPT screen output at a default level
 *
 * @param[in] p true or false
 */
void ipopt::screen_output(const bool p) {m_screen_out = p;}


std::string ipopt::human_readable_extra() const
{
	std::ostringstream s;
	s << "IPOPT - Max Iterations: " << m_max_iter << ", Tolerance: "<<m_tol<< ", Objective Function Tolerance: "<<m_acceptable_obj_change_tol << std::endl;
	return s.str();
}

}} //namespaces
