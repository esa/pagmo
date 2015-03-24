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

#ifndef IPOPT_PROBLEM_H
#define IPOPT_PROBLEM_H

#include <coin/IpTNLP.hpp>
#include "../../population.h"
#include "../../types.h"
#include "boost/array.hpp"


//Interface between Ipopt NLP and PaGMO problem

class ipopt_problem : public Ipopt::TNLP
{
public:
	/** default constructor */
	ipopt_problem(pagmo::population *);

	/** default destructor */
	virtual ~ipopt_problem();

	/**@name Overloaded from Ipopt::TNLP */
	//@{
	/** Method to return some info about the nlp */
	virtual bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
				  Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style);

	/** Method to return the bounds for my problem */
	virtual bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
				     Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u);

	/** Method to return the starting point for the algorithm */
	virtual bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
					bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
					Ipopt::Index m, bool init_lambda,
					Ipopt::Number* lambda);

	/** Method to return the objective value */
	virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value);

	/** Method to return the gradient of the objective */
	virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f);

	/** Method to return the constraint residuals */
	virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g);

	/** Method to return:
	 *   1) The structure of the jacobian (if "values" is NULL)
	 *   2) The values of the jacobian (if "values" is not NULL)
	 */
	virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
				Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
				Ipopt::Number* values);

	//@}

	/** @name Solution Methods */
	//@{
	/** This method is called when the algorithm is complete so the TNLP can store/write the solution */
	virtual void finalize_solution(Ipopt::SolverReturn status,
				       Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U,
				       Ipopt::Index m, const Ipopt::Number* g, const Ipopt::Number* lambda,
				       Ipopt::Number obj_value,
				       const Ipopt::IpoptData* ip_data,
				       Ipopt::IpoptCalculatedQuantities* ip_cq);
	//@}

private:
	/** @name Methods to block default compiler methods.
	 */
	//@{
	ipopt_problem();
	ipopt_problem(const ipopt_problem&);
	ipopt_problem& operator=(const ipopt_problem&);
	//@}
	//Points to a PaGMO population
	::pagmo::population *m_pop;
	//Number of non-zero entries in the Jacbian
	::Ipopt::Index len_jac;
	//Sparse representation of the Jacobian
	std::vector< ::Ipopt::Index> iJfun,jJvar;
	//Contains the variables that effect the objective function
	std::vector< ::Ipopt::Index> affects_obj;
	//Sorting criteria for the iJfun, jJvar entries to achieve constraint cache efficiency
	static bool cache_efficiency_criterion(boost::array<int,2>,boost::array<int,2>);
	// Internal caches used during evolution.
	::pagmo::decision_vector dv;
	::pagmo::fitness_vector fit;
	::pagmo::constraint_vector con;
};


#endif
