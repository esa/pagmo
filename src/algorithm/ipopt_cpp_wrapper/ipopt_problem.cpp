// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: MyNLP.cpp 1241 2008-05-29 23:05:26Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

#include "ipopt_problem.h"
#include <algorithm>

using namespace Ipopt;




/* Constructor. */
ipopt_problem::ipopt_problem(const ::pagmo::population &pop) : m_pop(pop)
{
	//We size the various members
	affects_obj.resize(0);
	dv.resize(m_pop.problem().get_dimension());
	fit.resize(m_pop.problem().get_f_dimension());
	con.resize(m_pop.problem().get_c_dimension());

	//We need to initialize the class members len_jac, iJfun, jJvar. We do this by first
	//storing the information in duples, ordering using the std::sort algorithm and then assigning
	//them to iJfun, jJvar. This is done to increase the constraint cache hits
	int lenG;
	std::vector<int> iGfun,jGvar;
	len_jac=0;
	std::vector< boost::array<int,2> > duples(0);
	boost::array<int,2> tmp;


	//If the problem has its set_sparsity implemented, store the values relevant to the constraints
	try {
		m_pop.problem().set_sparsity(lenG,iGfun,jGvar);
		for (int i = 0; i<lenG; ++i)
		{
			if (iGfun[i]!=0) //objective function gradient sparsity is not used by ipopt
			{
				tmp[0] = iGfun[i];
				tmp[1] = jGvar[i];
				duples.push_back(tmp);
				len_jac++;
			}
			else
			{
				affects_obj.push_back(jGvar[i]);
			}
		}
		//We now reorder the entries so that the cache will be hit avoiding useless
		//re-evaluations of the constraints in eval_jac_g when performing finite differences
		std::sort (duples.begin(), duples.end(), cache_efficiency_criterion);
	}
	//Otherwise assume no sparsity
	catch (not_implemented_error)
	{
		for (pagmo::problem::base::size_type j=0;j<m_pop.problem().get_dimension();++j)
		{
			for (pagmo::problem::base::size_type i=0;i<m_pop.problem().get_c_dimension();++i)
			{
				tmp[0] = i;
				tmp[1] = j;
				duples.push_back(tmp);
				len_jac++;
			}
		}
	}

	//And we finally assign the ordered duples to iJfun, jJvar
	iJfun.resize(len_jac);
	jJvar.resize(len_jac);
	for (int i = 0;i<len_jac;++i)
	{
		iJfun[i] = duples[i][0];
		jJvar[i] = duples[i][1];
	}
}

ipopt_problem::~ipopt_problem()
{}

//This function is the sort criteria for the elements of the sparse matrix J. First the variables,
//then the constraint number. i.e.
//Before sorting:
//iJfun = [0,3,0,0,2,1]
//jJvar = [1,2,0,2,1,0]
//After sorting:
//iJfun = [0,1,0,2,0,3]
//jJvar = [0,0,1,1,2,2]
bool ipopt_problem::cache_efficiency_criterion(boost::array<int,2> one,boost::array<int,2> two)
{
	if (one[1] < two[1]) {
		return true;
	}
	else {
		if (one[1] > two[1]) return false;
		else{ //they are equal!!! look to the other element
			if (one[0] < two[0]) return true;
			else return false;
		}
	}
}

bool ipopt_problem::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
				 Index& nnz_h_lag, IndexStyleEnum& index_style)
{
	// The problem dimension goes here
	n = m_pop.problem().get_dimension();

	// The number of constraints goes here.
	m = m_pop.problem().get_c_dimension();

	// Number of non zero entries in the jacobian
	nnz_jac_g = len_jac;

	// We use the standard c index style for row/col entries
	index_style = C_STYLE;

	return true;
}

bool ipopt_problem::get_bounds_info(Index n, Number* x_l, Number* x_u,
				    Index m, Number* g_l, Number* g_u)
{
	//Bounds on the decision vector
	for (pagmo::problem::base::size_type i=0; i<n;++i)
	{
		x_l[i] = m_pop.problem().get_lb()[i];
		x_u[i] = m_pop.problem().get_ub()[i];
	}

	//Bounds on equality constraints
	for (pagmo::problem::base::size_type i=0; i < (m-m_pop.problem().get_ic_dimension());++i)
	{
		g_l[i] = g_u[i] = 0.0;
	}

	//Bounds on inequality constraints
	for (pagmo::problem::base::size_type i = ( m - m_pop.problem().get_ic_dimension()); i < m;++i)
	{
		g_l[i] = - std::numeric_limits<double>::max();
		g_u[i] = 0.0;
	}



	return true;
}

bool ipopt_problem::get_starting_point(Index n, bool init_x, Number* x,
				       bool init_z, Number* z_L, Number* z_U,
				       Index m, bool init_lambda,
				       Number* lambda)
{
	// Here, we assume we only have starting values for x, if you code
	// your own NLP, you can provide starting values for the others if
	// you wish.
	assert(init_x == true);
	assert(init_z == false);
	assert(init_lambda == false);

	// we initialize x as the best individual of the population
	pagmo::population::size_type bestidx = m_pop.get_best_idx();
	std::copy(m_pop.get_individual(bestidx).cur_x.begin(),m_pop.get_individual(bestidx).cur_x.end(),x);

	return true;
}

bool ipopt_problem::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
	// return the value of the objective function
	std::copy(x,x+n,dv.begin());
	m_pop.problem().objfun(fit,dv);
	obj_value = fit[0];
	return true;
}

bool ipopt_problem::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
	double central_diff;
	const double h = 1e-8;
	std::copy(x,x+n,dv.begin());
	for (pagmo::decision_vector::size_type i=0; i<dv.size();++i)
	{
		grad_f[i] = 0;
	}

	for (size_t i =0;i<affects_obj.size();++i)
	{
		dv[affects_obj[i]] = dv[affects_obj[i]] + h;
		m_pop.problem().objfun(fit,dv);
		central_diff = fit[0];
		dv[affects_obj[i]] = dv[affects_obj[i]] - 2*h;
		m_pop.problem().objfun(fit,dv);
		central_diff = (central_diff-fit[0]) / 2 / h;
		grad_f[affects_obj[i]] = central_diff;
	}

	return true;
}

bool ipopt_problem::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
	// return the value of the constraints: g(x)
	std::copy(x,x+n,dv.begin());
	m_pop.problem().compute_constraints(con,dv);
	std::copy(con.begin(),con.end(),g);
	return true;
}

bool ipopt_problem::eval_jac_g(Index n, const Number* x, bool new_x,
			       Index m, Index nele_jac, Index* iRow, Index *jCol,
			       Number* values)
{
	if (values == NULL) {
		// return the structure of the jacobian of the constraints
		for (int i=0;i<len_jac;++i)
		{
			iRow[i] = iJfun[i];
			jCol[i] = jJvar[i];
		}
	}
	else {
		double central_diff;
		const double h = 1e-8;
		double mem;
		std::copy(x,x+n,dv.begin());
		for (int i=0;i<len_jac;++i)
		{
			mem = dv[jJvar[i]];
			dv[jJvar[i]] = dv[jJvar[i]] + h;
			m_pop.problem().compute_constraints(con,dv);
			central_diff = con[iJfun[i]];
			dv[jJvar[i]] = dv[jJvar[i]] - 2 * h;
			m_pop.problem().compute_constraints(con,dv);
			central_diff = (central_diff - con[iJfun[i]])/2/h;
			values[i] = central_diff;
			dv[jJvar[i]] = mem;
		}
	}

	return true;
}


void ipopt_problem::finalize_solution(SolverReturn status,
				      Index n, const Number* x, const Number* z_L, const Number* z_U,
				      Index m, const Number* g, const Number* lambda,
				      Number obj_value,
				      const IpoptData* ip_data,
				      IpoptCalculatedQuantities* ip_cq)
{
	// here is where we would store the solution to variables, or write to a file, etc
	// so we could use the solution.
}
