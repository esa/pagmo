/*****************************************************************************
  *   Copyright (C) 2004-2015 The PaGMO development team,                     *
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
 * Allows to specify some of the parameters of the IPOPT solver. The algorithm stopping criteria will
 * be 1) if the number of iteration exceeds max_iter 2) the three tolerances are met
 *
 * @param[in] max_iter Maximum number of major iterations (refer to IPOPT manual).
 * @param[in] constr_viol_tol Constraint violation tolerance
 * @param[in] dual_inf_tol Dual infeasibility tolerance
 * @param[in] compl_inf_tol Complementary feasibility tolerance
 * @param[in] nlp_scaling_method Select if the "gradient-based" scaling of the  NLP should be used
 * @param[in] obj_scaling_factor Scaling factor for the objective function.
 * @param[in] mu_init Initial value for the barrier parameter.
 * @throws value_error if max_iter or tolerances are negative
 */

ipopt::ipopt(	const int &max_iter,
		const double &constr_viol_tol, 
		const double &dual_inf_tol, 
		const double &compl_inf_tol,
		const bool &nlp_scaling_method,
		const double &obj_scaling_factor,
		const double &mu_init) :
		m_max_iter(max_iter),m_constr_viol_tol(constr_viol_tol),
		m_dual_inf_tol(dual_inf_tol), m_compl_inf_tol(compl_inf_tol), 
		m_nlp_scaling_method(nlp_scaling_method), m_obj_scaling_factor(obj_scaling_factor), m_mu_init(mu_init)
		
{
	if (max_iter < 0) {
		pagmo_throw(value_error,"number of maximum iterations cannot be negative");
	}
	if (constr_viol_tol < 0) {
		pagmo_throw(value_error,"tolerance is not in ]0,1[");
	}
	if (dual_inf_tol < 0) {
		pagmo_throw(value_error,"obj_tol is not in ]0,1[");
	}
	if (compl_inf_tol < 0) {
		pagmo_throw(value_error,"obj_tol is not in ]0,1[");
	}
	if (mu_init <= 0) {
		pagmo_throw(value_error,"mu_init is not in ]0, inf[");
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
 * At the end, the velocity is updated
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
	::Ipopt::SmartPtr< ::Ipopt::IpoptApplication> m_app = new ::Ipopt::IpoptApplication(m_screen_output,false);

	if (!m_nlp_scaling_method) {
		m_app->Options()->SetStringValue("nlp_scaling_method", "none");
	}
	m_app->Options()->SetStringValue("hessian_approximation", "limited-memory");
	m_app->Options()->SetIntegerValue("print_level", 5);

	// Termination Criteria 1: iterations
	m_app->Options()->SetIntegerValue("max_iter", m_max_iter);

	// Termination Criteria 2: tolerance
	m_app->Options()->SetNumericValue("tol", 1.);
	m_app->Options()->SetNumericValue("dual_inf_tol", m_dual_inf_tol);
	m_app->Options()->SetNumericValue("constr_viol_tol", m_constr_viol_tol);
	m_app->Options()->SetNumericValue("compl_inf_tol", m_compl_inf_tol);
	m_app->Options()->SetNumericValue("obj_scaling_factor", m_obj_scaling_factor);
	m_app->Options()->SetNumericValue("mu_init", m_mu_init);



	// Intialize the IpoptApplication and process the options
	Ipopt::ApplicationReturnStatus status;
	status = m_app->Initialize();
	if (status != Ipopt::Solve_Succeeded) {
		pagmo_throw(value_error, "Error during IPOPT initialization!");
	}

	// Ask Ipopt to solve the problem
	status = m_app->OptimizeTNLP(pagmo_nlp);
}

/// Algorithm name
std::string ipopt::get_name() const
{
	return "IPOPT";
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string ipopt::human_readable_extra() const
{
	std::ostringstream s;
	s << "major_iter:" << m_max_iter << " ";
	s << "constr_viol_tol:"<< m_constr_viol_tol<<" ";
	s << "dual_inf_tol:"<< m_dual_inf_tol<<" ";
	s << "compl_inf_tol:"<< m_compl_inf_tol<<" ";
	s << "nlp_scaling_method:"<< (m_nlp_scaling_method? "True":"False") << " ";
	s << "obj_scaling_factor:"<< m_obj_scaling_factor<<" ";
	s << "mu_init:"<< m_mu_init;
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::ipopt)
