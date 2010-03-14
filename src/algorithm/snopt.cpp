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
#include "snopt.h"
#include "snopt_cpp_wrapper/snopt_PAGMO.h"
#include "snopt_cpp_wrapper/snfilewrapper_PAGMO.h"
#include "snopt_cpp_wrapper/snoptProblem_PAGMO.h"

//These are here to solve a possible problem with the shared f2c libraries during linking
extern "C"{
	void MAIN_(){}
	void MAIN__(){}
	void _MAIN_(){}
}


///This is a static function that is used to wrap the snopt libraries. It has the propotype that is required
///and uses the memory location pointed by cu to read all informations about the PaGMO problem.
static int snopt_function_(integer    *Status, integer *n,    doublereal *x,
			   integer    *needF,  integer *neF,  doublereal *F,
			   integer    *needG,  integer *neG,  doublereal *G,
			   char       *cu,     integer *lencu,
			   integer    *iu,    integer *leniu,
			   doublereal *ru,    integer *lenru )
{
	//1 - We retrieve the pointer to the base problem (PaGMO) we have 'hidden' in *cu
	pagmo::problem::base *prob;
	prob = (pagmo::problem::base*)cu;

	//2 - We retrieve the pointer to the allocated memory area containing the std::vector where
	//to copy the snopt x[] array
	pagmo::decision_vector* chromosome = (pagmo::decision_vector*)ru;
	for (size_t i = 0;i < ( prob->get_dimension() - prob->get_i_dimension() );i++) (*chromosome)[i] = x[i];

	//We finally assign to F[0] (the objective function) the value returned by the problem
	pagmo::fitness_vector fit(1);
	try {
		prob->objfun(fit,*chromosome);
		F[0] = fit[0];
		}
	catch (value_error) {
		*Status = -1; //signals to snopt that the evaluation of the objective function had numerical difficulties
	}
	//And to F[.] the constraint values (equality first)
	std::vector<double> con(prob->get_c_dimension());
	prob->compute_constraints(con, *chromosome);
	for (int i=0;i<prob->get_c_dimension();++i) F[i+1] = con[i];

	return 0;
}

namespace pagmo { namespace algorithm {
/// Constructor
snopt::snopt(const int major,const double feas, const double opt) : m_major(major),m_feas(feas),m_opt(opt)
{
	if (major < 0) {
		pagmo_throw(value_error,"number of major iterations cannot be negative");
	}
	if (feas < 0 || feas > 1) {
		pagmo_throw(value_error,"feasibility cireria ");
	}
	if (opt < 0 || opt > 1) {
		pagmo_throw(value_error,"number of major iterations cannot be negative");
	}

}
/// Clone method.
base_ptr snopt::clone() const
{
	return base_ptr(new snopt(*this));
}

/// Evolve implementation.
/**
 * Run SNOPT for the number of iterations specified in the constructors.
 * At the end velocity is updated
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void snopt::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type D = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension(), prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = pop.size();
	const problem::base::size_type Dc = D - prob_i_dimension;
	const std::vector<double>::size_type D_ineqc = prob.get_ic_dimension();
	const std::vector<double>::size_type D_eqc = prob_c_dimension - D_ineqc;

	//We perform some checks to determine wether the problem/population are suitable for DE
	if ( prob_i_dimension != 0  ) {
		pagmo_throw(value_error,"No integer part allowed yet....");
	}

	if ( Dc == 0  ) {
		pagmo_throw(value_error,"No continuous part....");
	}

	if ( prob_f_dimension != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and DE is not suitable to solve it");
	}

	// Get out if there is nothing to do.
	if (NP == 0 || m_major == 0) {
		return;
	}

	// We allocate memory for the decision vector that will be used in the snopt_function_
	di_comodo.resize(Dc);

	// We construct a SnoptProblem_PAGMO passing the pointers to the problem and the allocated
	//memory area for the di_comodo vector
	snoptProblem_PAGMO SnoptProblem(prob, &di_comodo);

	// Allocate and initialize;
	integer n     =  Dc;

	// Box-constrained non-linear optimization
	integer neF   =  1 + prob_c_dimension;
	integer lenA  = 1; //needs to be at least 1

	//No linear part lenA = 0
	integer *iAfun = new integer[lenA];
	integer *jAvar = new integer[lenA];
	doublereal *A  = new doublereal[lenA];

	//Pattern structure as defined by the problem
	integer lenG   = prob.get_neG();
	integer *iGfun = new integer[lenG];
	integer *jGvar = new integer[lenG];
	for (int i=0;i<lenG;i++){
		iGfun[i] = prob.get_iGfun()[i];
		jGvar[i] = prob.get_jGvar()[i];
	}

	//Decision vector memory allocation
	doublereal *x      = new doublereal[n];
	doublereal *xlow   = new doublereal[n];
	doublereal *xupp   = new doublereal[n];
	doublereal *xmul   = new doublereal[n];
	integer    *xstate = new    integer[n];

	//Objective function memory allocation
	doublereal *F      = new doublereal[neF];
	doublereal *Flow   = new doublereal[neF];
	doublereal *Fupp   = new doublereal[neF];
	doublereal *Fmul   = new doublereal[neF];
	integer    *Fstate = new integer[neF];

	integer nxnames = 1;
	integer nFnames = 1;
	char *xnames = new char[nxnames*8];
	char *Fnames = new char[nFnames*8];

	integer    ObjRow = 0;
	doublereal ObjAdd = 0;

	// Set the upper and lower bounds. And The initial Guess
	int bestidx = pop.get_best_idx();
	for (int i = 0; i < n; i++){
		xlow[i]   = lb[i];
		xupp[i]   = ub[i];
		xstate[i] =    0;
		x[i] = 1;//pop.get_individual(bestidx).cur_x[i];
	}

	// Set the bounds for objective, equality and inequality constraints
	Flow[0] = -1e20;
	Fupp[0] = 1e20;
	F[0] = pop.get_individual(bestidx).cur_f[0];
	for (int i=0;i<D_eqc;++i) {
		Flow[i+1] = 0;
		Fupp[i+1] = 0;
	}
	for (int i=0;i<D_ineqc;++i) {
		Flow[i+1+D_eqc] = -1e20;
		Fupp[i+1+D_eqc] = 0;
	}

	// Load the data for SnoptProblem ...
	char *tmp = "Toy0.out";
	SnoptProblem.setPrintFile  ( tmp );
	SnoptProblem.setProblemSize( n, neF );
	SnoptProblem.setObjective  ( ObjRow, ObjAdd );
	SnoptProblem.setA          ( lenA, iAfun, jAvar, A );
	SnoptProblem.setG          ( lenG, iGfun, jGvar );
	SnoptProblem.setX          ( x, xlow, xupp, xmul, xstate );
	SnoptProblem.setF          ( F, Flow, Fupp, Fmul, Fstate );
	SnoptProblem.setXNames     ( xnames, nxnames );
	SnoptProblem.setFNames     ( Fnames, nFnames );
	SnoptProblem.setProbName   ( "Toy0" );
	SnoptProblem.setUserFun    ( snopt_function_ );

	// snopta will compute the Jacobian by finite-differences.
	// The user has the option of calling  snJac  to define the
	// coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
	SnoptProblem.computeJac    ();
	SnoptProblem.setIntParameter( "Derivative option", 0 );
	SnoptProblem.setIntParameter( "Major iterations limit", m_major);
	SnoptProblem.setRealParameter( "Major feasibility tolerance", m_feas);
	SnoptProblem.setRealParameter( "Major optimality tolerance", m_opt);

	integer Cold = 0, Basis = 1, Warm = 2;

	//HERE WE CALL snoptA routine!!!!!
	SnoptProblem.solve( Cold );

	//Save the final point
	for (integer i=0;i<n;i++) di_comodo[i] = x[i];
	decision_vector newx = di_comodo;
	std::transform(di_comodo.begin(), di_comodo.end(), pop.get_individual(bestidx).cur_x.begin(), di_comodo.begin(),std::minus<double>());
	pop.set_x(bestidx,newx);
	pop.set_v(bestidx,di_comodo);

	//Clean up memory allocated to call the snoptA routine
	delete []iAfun;  delete []jAvar;  delete []A;
	delete []iGfun;  delete []jGvar;

	delete []x;      delete []xlow;   delete []xupp;
	delete []xmul;   delete []xstate;

	delete []F;      delete []Flow;   delete []Fupp;
	delete []Fmul;   delete []Fstate;

	delete []xnames; delete []Fnames;

}

std::string snopt::human_readable_extra() const
{
	std::ostringstream s;
	s << "SNOPT - Major Iterations: " << m_major << ", Feasibility: "<<m_feas<< ", Optimality: "<<m_opt << std::endl;
	return s.str();
}

}} //namespaces
