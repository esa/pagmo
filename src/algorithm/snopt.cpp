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
	decision_vector* chromosome = (decision_vector*)ru;
	for (size_t i = 0;i < ( prob->get_dimension() - prob->get_i_dimension() );i++) (*chromosome)[i] = x[i];

	//We finally assign to F[0] (the objective function) the value returned by the problem
	fitness_vector fit;
	try {
		prob->objfun(fit,chromosome);
		F[0] = fit[0];
		}
	catch (value_error) {
		*Status = -1; //signals to snopt that the evaluation of the objective function had numerical difficulties
	}
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

void snopt::evolve(const population &pop) const
{
	// We allocate memory for the decision vector that will be used in the snopt_function_
	di_comodo.resize(pop.problem().getDimension());

	// We construct a SnoptProblem_PAGMO passing the pointers to the problem and the allocated
	//memory area for the di_comodo vector
	snoptProblem_PAGMO SnoptProblem(pop.problem(), &di_comodo);

	// Allocate and initialize;
	integer n     =  pop.problem().getDimension();

	// Box-constrained non-linear optimization
	integer neF   =  1;
	integer lenA  = 0;

	//No linear part lenA = 0
	integer *iAfun = new integer[lenA];
	integer *jAvar = new integer[lenA];
	doublereal *A  = new doublereal[lenA];

	//Pattern structure as defined in the GOproblem
	integer lenG   = pop.problem().get_neG();
	integer *iGfun = new integer[lenG];
	integer *jGvar = new integer[lenG];
	for (int i=0;i<lenG;i++){
		iGfun[i] = pop.problem().get_iGfun()[i];
		jGvar[i] = pop.problem().get_jGvar()[i];
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
	for (int i = 0; i < n; i++){
		xlow[i]   =  pop.problem().getLB()[i];
		xupp[i]   = pop.problem().getUB()[i];
		xstate[i] =    0;
		x[i] = pop.extractBestIndividual().getDecisionVector()[i];
	}

	Flow[0] = -1e20; //Flow[1] = -1e20; Flow[2] = -1e20;
	Fupp[0] =  1e20;// Fupp[1] =   4.0; Fupp[2] =  5.0;
	F[0] = pop.extractBestIndividual().getFitness();


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
	SnoptProblem.setIntParameter( "Major iterations limit", major_iterations);
	SnoptProblem.setRealParameter( "Major feasibility tolerance", feasibility);
	SnoptProblem.setRealParameter( "Major optimality tolerance", optimality);

	integer Cold = 0, Basis = 1, Warm = 2;

	//HERE WE CALL snoptA routine!!!!!
	SnoptProblem.solve( Cold );

	//Save the final point
	std::vector<double> retval(n);
	for (integer i=0;i<n;i++) retval[i] = x[i];
	double fitness = F[0];

	//Save the velocity (unchanged) and the old fitness
	std::vector<double> vel = pop.extractBestIndividual().getVelocity();
	double oldfit =  pop.extractBestIndividual().getFitness();

	//Clean up memory allocated to call the snoptA routine
	delete []iAfun;  delete []jAvar;  delete []A;
	delete []iGfun;  delete []jGvar;

	delete []x;      delete []xlow;   delete []xupp;
	delete []xmul;   delete []xstate;

	delete []F;      delete []Flow;   delete []Fupp;
	delete []Fmul;   delete []Fstate;

	delete []xnames; delete []Fnames;

	//Construct the new population (same as old except best_individual if improved)
	Population NewPop(pop);
	if (fitness < oldfit) NewPop.replace_best(Individual(retval,vel,fitness));
	return (NewPop);
}

void snopt_algorithm::log(std::ostream &s) const
{
	s << "SNOPT - Major Iterations: " << major_iterations << ", Feasibility: "<<feasibility<< ", Optimality: "<<optimality << std::endl;
}

}} //namespaces
