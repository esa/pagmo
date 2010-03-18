// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: MyNLP.cpp 1241 2008-05-29 23:05:26Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

#include "ipopt_problem.h"
#include "../../problem/base.h"

using namespace Ipopt;

/* Constructor. */
ipopt_problem::ipopt_problem(const ::pagmo::problem::base prob) : m_prob(prob*)
{
	//We need here to initialize the class members len_jac, iJfun, jJvar. We do this by first
	//storing the information in duples, ordering it using the std::sort algorithm and then assigning
	//them to iJfun, jJvar. This is done to increase the constraint cache hits
	int lenG;
	std::vector<double> iGfun,jGvar;
	len_jac=0;
	std::vector< boost::array<int,2> > duples(0);
	boost::array<int,2> tmp;
	//If the problem has its set_sparsity implemented, store the values relevant to the constraints
	try {
		m_prob->set_sparsity(lenG,iGfun,jGvar);
		for (int i = 0; i<lenG; ++i)
		{
			if (iGfun!=0) //objective function gradient sparsity is not used by ipopt
			{
				tmp[0] = iGfun[i];
				tmp[1] = jGvar[i];
				duples.push_back(tmp);
				len_jac++;
			}
		}
		//We now reorder the entries so that the cache will be hit avoiding useless
		//re-evaluations of the constraints in eval_jac_g when performing finite differences
		std::sort (duples.begin(), duples.end(), cache_efficiency_criterion);
	}
	//Otherwise assume no sparsity
	catch (not_implemented_error)
	{
		for (pagmo::problem::base::size_type j=0;j<m_prob->get_dimension();++j)
		{
			for (pagmo::problem::base::size_type i=0;i<m_prob->get_c_dimension();++i)
			{
				tmp[0] = i;
				tmp[1] = j;
				duples.push_back(tmp);
				len_jac++;
			}
		}
	}

	//And we finally assign the ordered duples to iJfun, jJvar
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
//[0,3,0,0,2,1]
//[1,2,0,2,1,0]
//After sorting:
//[0,1,0,2,0,3]
//[0,0,1,1,2,2]
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
	n = m_prob->get_dimension();

	// The number of constraints goes here. Each equality constraint is transformed into two inequalities
	m = 2 * (m_prob->get_c_dimension() - m_prob->get_ic_dimension()) + m_prob->get_ic_dimension();

	// Number of non zero entries in the jacobian
	nnz_jac_g = len_jac;

	// We use the standard c index style for row/col entries
	index_style = C_STYLE;

	return true;
}

bool ipopt_problem::get_bounds_info(Index n, Number* x_l, Number* x_u,
				    Index m, Number* g_l, Number* g_u)
{
	// here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
	// If desired, we could assert to make sure they are what we think they are.
	assert(n == 2);
	assert(m == 1);

	// x1 has a lower bound of -1 and an upper bound of 1
	x_l[0] = -1.0;
	x_u[0] = 1.0;

	// x2 has no upper or lower bound, so we set them to
	// a large negative and a large positive number.
	// The value that is interpreted as -/+infinity can be
	// set in the options, but it defaults to -/+1e19
	x_l[1] = -1.0e19;
	x_u[1] = +1.0e19;

	// we have one equality constraint, so we set the bounds on this constraint
	// to be equal (and zero).
	g_l[0] = g_u[0] = 0.0;

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

	// we initialize x in bounds, in the upper right quadrant
	x[0] = 0.5;
	x[1] = 1.5;

	return true;
}

bool ipopt_problem::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
	// return the value of the objective function
	Number x2 = x[1];
	obj_value = -(x2 - 2.0) * (x2 - 2.0);
	return true;
}

bool ipopt_problem::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
	// return the gradient of the objective function grad_{x} f(x)

	// grad_{x1} f(x): x1 is not in the objective
	grad_f[0] = 0.0;

	// grad_{x2} f(x):
	Number x2 = x[1];
	grad_f[1] = -2.0*(x2 - 2.0);

	return true;
}

bool ipopt_problem::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
	// return the value of the constraints: g(x)
	Number x1 = x[0];
	Number x2 = x[1];

	g[0] = -(x1*x1 + x2 - 1.0);

	return true;
}

bool ipopt_problem::eval_jac_g(Index n, const Number* x, bool new_x,
			       Index m, Index nele_jac, Index* iRow, Index *jCol,
			       Number* values)
{
	if (values == NULL) {
		// return the structure of the jacobian of the constraints

		// element at 1,1: grad_{x1} g_{1}(x)
		iRow[0] = 1;
		jCol[0] = 1;

		// element at 1,2: grad_{x2} g_{1}(x)
		iRow[1] = 1;
		jCol[1] = 2;
	}
	else {
		// return the values of the jacobian of the constraints
		Number x1 = x[0];

		// element at 1,1: grad_{x1} g_{1}(x)
		values[0] = -2.0 * x1;

		// element at 1,2: grad_{x1} g_{1}(x)
		values[1] = -1.0;
	}

	return true;
}

bool ipopt_problem::eval_h(Index n, const Number* x, bool new_x,
			   Number obj_factor, Index m, const Number* lambda,
			   bool new_lambda, Index nele_hess, Index* iRow,
			   Index* jCol, Number* values)
{
	if (values == NULL) {
		// return the structure. This is a symmetric matrix, fill the lower left
		// triangle only.

		// element at 1,1: grad^2_{x1,x1} L(x,lambda)
		iRow[0] = 1;
		jCol[0] = 1;

		// element at 2,2: grad^2_{x2,x2} L(x,lambda)
		iRow[1] = 2;
		jCol[1] = 2;

		// Note: off-diagonal elements are zero for this problem
	}
	else {
		// return the values

		// element at 1,1: grad^2_{x1,x1} L(x,lambda)
		values[0] = -2.0 * lambda[0];

		// element at 2,2: grad^2_{x2,x2} L(x,lambda)
		values[1] = -2.0 * obj_factor;

		// Note: off-diagonal elements are zero for this problem
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
	// so we could use the solution. Since the solution is displayed to the console,
	// we currently do nothing here.
}
