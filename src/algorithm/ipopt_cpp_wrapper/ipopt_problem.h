// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: MyNLP.hpp 949 2007-03-27 00:41:26Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

#ifndef IPOPT_PROBLEM_H
#define IPOPT_PROBLEM_H

#include <coin/IpTNLP.hpp>
#include "../../population.cpp"
#include "../../types.h"
#include "boost/array.hpp"

using namespace Ipopt;

//Interface between Ipopt NLP and PaGMO problem

class ipopt_problem : public TNLP
{
public:
	/** default constructor */
	ipopt_problem(const pagmo::population &);

	/** default destructor */
	virtual ~ipopt_problem();

	/**@name Overloaded from Ipopt::TNLP */
	//@{
	/** Method to return some info about the nlp */
	virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
				  Index& nnz_h_lag, IndexStyleEnum& index_style);

	/** Method to return the bounds for my problem */
	virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
				     Index m, Number* g_l, Number* g_u);

	/** Method to return the starting point for the algorithm */
	virtual bool get_starting_point(Index n, bool init_x, Number* x,
					bool init_z, Number* z_L, Number* z_U,
					Index m, bool init_lambda,
					Number* lambda);

	/** Method to return the objective value */
	virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

	/** Method to return the gradient of the objective */
	virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

	/** Method to return the constraint residuals */
	virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

	/** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
	virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
				Index m, Index nele_jac, Index* iRow, Index *jCol,
				Number* values);

	//@}

	/** @name Solution Methods */
	//@{
	/** This method is called when the algorithm is complete so the TNLP can store/write the solution */
	virtual void finalize_solution(SolverReturn status,
				       Index n, const Number* x, const Number* z_L, const Number* z_U,
				       Index m, const Number* g, const Number* lambda,
				       Number obj_value,
				       const IpoptData* ip_data,
				       IpoptCalculatedQuantities* ip_cq);
	//@}

private:
	/**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually 
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *  
   */
	//@{
	ipopt_problem();
	ipopt_problem(const ipopt_problem&);
	ipopt_problem& operator=(const ipopt_problem&);
	//@}
	//Points to a PaGMO population
	const ::pagmo::population m_pop;
	//Number of non-zero entries in the Jacbian
	int len_jac;
	//Sparse representation of the Jacobian
	std::vector<int> iJfun,jJvar;
	//Contains the variables that effect the objective function
	std::vector<int> affects_obj;
	// Internal caches used during evolution.
	::pagmo::decision_vector dv;
	::pagmo::fitness_vector fit;
	::pagmo::constraint_vector con;
};


#endif
