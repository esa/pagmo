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


//This is a static function that is used to wrap the snopt libraries. It has the propotype that is required
//and uses the memory location pointed by cu to read all informations about the PaGMO problem.
static int snopt_function_(integer    *Status, integer *n,    doublereal *x,
			   integer    *needF,  integer *neF,  doublereal *F,
			   integer    *needG,  integer *neG,  doublereal *G,
			   char       *cu,     integer *lencu,
			   integer    *iu,    integer *leniu,
			   doublereal *ru,    integer *lenru )
{
	(void)n;
	(void)needF;
	(void)neF;
	(void)needG;
	(void)neG;
	(void)G;
	(void)lencu;
	(void)iu;
	(void)leniu;
	(void)lenru;
	//1 - We retrieve the pointer to the base problem (PaGMO) we have 'hidden' in *cu
	pagmo::problem::base *prob;
	prob = (pagmo::problem::base*)cu;

	//2 - We retrieve the pointer to the allocated memory area containing the std::vector where
	//to copy the snopt x[] array
	pagmo::algorithm::snopt::preallocated_memory* preallocated = (pagmo::algorithm::snopt::preallocated_memory*)ru;
	for (size_t i = 0;i < ( prob->get_dimension() - prob->get_i_dimension() );i++) (preallocated->x)[i] = x[i];

	//1 - to F[0] (the objective function) the value returned by the problem
	try {
		prob->objfun(preallocated->f,preallocated->x);
		F[0] = preallocated->f[0];
		}
	catch (value_error) {
		*Status = -1; //signals to snopt that the evaluation of the objective function had numerical difficulties
	}
	//2 - and to F[.] the constraint values (equalities first)
	try{
		prob->compute_constraints(preallocated->c, preallocated->x);
		for (pagmo::problem::base::size_type i=0;i<prob->get_c_dimension();++i) F[i+1] = preallocated->c[i];
	}
	catch (value_error) {
		*Status = -1; //signals to snopt that the evaluation of the objective function had numerical difficulties
	}

	return 0;
}

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify some of the parameters of the SNOPT solver.
 *
 * @param[in] major Number of major iterations (refer to SNOPT manual).
 * @param[in] feas Feasibility tolerance (refer to SNOPT manual).
 * @param[in] opt Optimality tolerance (refer to SNOPT manual).
 * @throws value_error if major is not positive, and feas,opt are not in \f$]0,1[\f$
 */
snopt::snopt(const int major,const double feas, const double opt) : m_major(major),m_feas(feas),m_opt(opt), m_file_out(false)
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
 * Run SNOPT with the parameters specified in the constructor
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
	const std::string name = prob.get_name();

	//We perform some checks to determine wether the problem/population are suitable for SNOPT
	if ( prob_i_dimension != 0  ) {
		pagmo_throw(value_error,"No integer part allowed yet....");
	}

	if ( Dc == 0  ) {
		pagmo_throw(value_error,"No continuous part....");
	}

	if ( prob_f_dimension != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and SNOPT is not suitable to solve it");
	}

	// Get out if there is nothing to do.
	if (NP == 0 || m_major == 0) {
		return;
	}

	// We allocate memory for the decision vector that will be used in the snopt_function_
	di_comodo.x.resize(Dc);
	di_comodo.c.resize(prob_c_dimension);
	di_comodo.f.resize(prob_f_dimension);


	// We construct a SnoptProblem_PAGMO passing the pointers to the problem and the allocated
	//memory area for the di_comodo vector
	snoptProblem_PAGMO SnoptProblem(prob, &di_comodo);

	// Allocate and initialize;
	integer n     =  Dc;

	// Box-constrained non-linear optimization
	integer neF   =  1 + prob_c_dimension;

	//Memory sizing of A
	integer lenA  = Dc * (1 + prob_c_dimension); //overestimate
	integer *iAfun = new integer[lenA];
	integer *jAvar = new integer[lenA];
	doublereal *A  = new doublereal[lenA];


	//Memory sizing of G
	int lenG =Dc * (1 + prob_c_dimension); //overestimate
	integer *iGfun = new integer[lenG];
	integer *jGvar = new integer[lenG];



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
	for (pagmo::problem::base::size_type i = 0; i < Dc; i++){
		xlow[i]   = lb[i];
		xupp[i]   = ub[i];
		xstate[i] =    0;
		x[i] = pop.get_individual(bestidx).cur_x[i];
	}

	// Set the bounds for objective, equality and inequality constraints
	// 1 - Objective function
	Flow[0] = -std::numeric_limits<double>::max();
	Fupp[0] = std::numeric_limits<double>::max();
	F[0] = pop.get_individual(bestidx).cur_f[0];
	// 2 - Equality constraints
	for (pagmo::problem::base::size_type i=0;i<D_eqc;++i) {
		Flow[i+1] = 0;
		Fupp[i+1] = 0;
	}
	// 3 - Inequality constraints
	for (pagmo::problem::base::size_type i=0;i<D_ineqc;++i) {
		Flow[i+1+D_eqc] = -std::numeric_limits<double>::max();
		Fupp[i+1+D_eqc] = 0;
	}

	// Load the data for SnoptProblem ...
	SnoptProblem.setProblemSize( n, neF );
	SnoptProblem.setNeG( lenG );
	SnoptProblem.setNeA( lenA );
	SnoptProblem.setA          ( lenA, iAfun, jAvar, A );
	SnoptProblem.setG          ( lenG, iGfun, jGvar );
	SnoptProblem.setObjective  ( ObjRow, ObjAdd );
	SnoptProblem.setX          ( x, xlow, xupp, xmul, xstate );
	SnoptProblem.setF          ( F, Flow, Fupp, Fmul, Fstate );
	SnoptProblem.setXNames     ( xnames, nxnames );
	SnoptProblem.setFNames     ( Fnames, nFnames );
	SnoptProblem.setProbName   ( name.c_str() ); //This is limited to be 8 characters!!!
	SnoptProblem.setUserFun    ( snopt_function_ );

	//We set some parameters
	if (m_screen_output) SnoptProblem.setIntParameter("Summary file",6);
	if (m_file_out)   SnoptProblem.setPrintFile   ( name.c_str() );
	SnoptProblem.setIntParameter ( "Derivative option", 0 );
	SnoptProblem.setIntParameter ( "Major iterations limit", m_major);
	SnoptProblem.setIntParameter ( "Iterations limit",100000);
	SnoptProblem.setRealParameter( "Major feasibility tolerance", m_feas);
	SnoptProblem.setRealParameter( "Major optimality tolerance", m_opt);


	//We set the sparsity structure
	int neG;
	try
	{
		std::vector<int> iGfun_vect, jGvar_vect;
		prob.set_sparsity(neG,iGfun_vect,jGvar_vect);
		for (int i=0;i < neG;i++)
		{
			iGfun[i] = iGfun_vect[i];
			jGvar[i] = jGvar_vect[i];
		}
		SnoptProblem.setNeG( neG );
		SnoptProblem.setNeA( 0 );
		SnoptProblem.setG( lenG, iGfun, jGvar );

	} //the user did implement the sparsity in the problem
	catch (not_implemented_error)
	{
		SnoptProblem.computeJac();
		neG = SnoptProblem.getNeG();
	} //the user did not implement the sparsity in the problem


	if (m_screen_output)
	{
		std::cout << "PaGMO 4 SNOPT:" << std::endl << std::endl;
		std::cout << "Sparsity pattern set, NeG = " << neG << std::endl;
		std::cout << "iGfun: [";
		for (int i=0; i<neG-1; ++i) std::cout << iGfun[i] << ",";
		std::cout << iGfun[neG-1] << "]" << std::endl;
		std::cout << "jGvar: [";
		for (int i=0; i<neG-1; ++i) std::cout << jGvar[i] << ",";
		std::cout << jGvar[neG-1] << "]" << std::endl;
	}

	integer Cold = 0;

	//HERE WE CALL snoptA routine!!!!!
	SnoptProblem.solve( Cold );

	//Save the final point making sure it is within the linear bounds
	std::copy(x,x+n,di_comodo.x.begin());
	decision_vector newx = di_comodo.x;
	std::transform(di_comodo.x.begin(), di_comodo.x.end(), pop.get_individual(bestidx).cur_x.begin(), di_comodo.x.begin(),std::minus<double>());
	for (integer i=0;i<n;i++)
	{
		newx[i] = std::min(std::max(lb[i],newx[i]),ub[i]);
	}

	pop.set_x(bestidx,newx);
	pop.set_v(bestidx,di_comodo.x);

	//Clean up memory allocated to call the snoptA routine
	delete []iAfun;  delete []jAvar;  delete []A;
	delete []iGfun;  delete []jGvar;

	delete []x;      delete []xlow;   delete []xupp;
	delete []xmul;   delete []xstate;

	delete []F;      delete []Flow;   delete []Fupp;
	delete []Fmul;   delete []Fstate;

	delete []xnames; delete []Fnames;

}


/// Activate file output
/**
 * Activate SNOPT screen output by setting iPrint to 15 and setting the file name to the problem typeid
 *
 * @param[in] p true or false
 */
void snopt::file_output(const bool p) {m_file_out = p;}

/// Algorithm name
std::string snopt::get_name() const
{
	return "SNOPT";
}


/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string snopt::human_readable_extra() const
{
	std::ostringstream s;
	s << "major_iter:" << m_major << " feas_tol:"<<m_feas<< " opt_tol:"<<m_opt << std::endl;
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::snopt)
