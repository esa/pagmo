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

#include <boost/numeric/conversion/cast.hpp>

#include "ipopt_problem.h"
#include "../../exceptions.h"

using namespace Ipopt;

/* Constructor. */
ipopt_problem::ipopt_problem(pagmo::population *pop) : m_pop(pop)
{
	//We size the various members
	affects_obj.resize(0);
	dv.resize(m_pop->problem().get_dimension());
	fit.resize(m_pop->problem().get_f_dimension());
	con.resize(m_pop->problem().get_c_dimension());

	//We need to initialize the class members len_jac, iJfun, jJvar. We do this by first
	//storing the information in duples, ordering using the std::sort algorithm and then assigning
	//them to iJfun, jJvar. This is done to increase the constraint cache hits
	::Ipopt::Index lenG;
	std::vector< ::Ipopt::Index> iGfun,jGvar;
	len_jac=0;
	std::vector< boost::array< ::Ipopt::Index,2> > duples(0);
	boost::array< ::Ipopt::Index,2> tmp;


	//If the problem has its set_sparsity implemented, store the values relevant to the constraints
	try {
		m_pop->problem().set_sparsity(lenG,iGfun,jGvar);
		for (::Ipopt::Index i = 0; i<lenG; ++i)
		{
			if (iGfun[i]!=0) //objective function gradient sparsity is not used by ipopt
			{
				tmp[0] = iGfun[i] - 1; //HERE WAS THE UNDEBUGGABLE!! -1 ffffffkkkkkk
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
		for (pagmo::problem::base::size_type j=0;j<m_pop->problem().get_dimension();++j)
		{
			affects_obj.push_back(j);
			for (pagmo::problem::base::size_type i=0;i<m_pop->problem().get_c_dimension();++i)
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
	for (::Ipopt::Index i = 0;i<len_jac;++i)
	{
		iJfun[i] = duples[i][0];
		jJvar[i] = duples[i][1];
		//std::cout << "[" << iJfun[i] << "," << jJvar[i] << "]" << std::endl;
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


bool ipopt_problem::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
				 Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style)
{
	(void)nnz_h_lag; //hessian is evaluated numerically using BFGS
	//Problem dimension
	n = m_pop->problem().get_dimension();

	//Dimension of the constraint vector (equality and inequality)
	m = m_pop->problem().get_c_dimension();

	// Number of non -zero elements in Jacobian
	nnz_jac_g = len_jac;

	// We use the standard fortran index style for row/col entries
	index_style = C_STYLE;

	return true;
}

bool ipopt_problem::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
					Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u)
{
	//Bounds on the decision vector
	for (::Ipopt::Index i=0; i<n;++i)
	{
		x_l[i] = m_pop->problem().get_lb()[i];
		x_u[i] = m_pop->problem().get_ub()[i];
	}

	//Bounds on equality constraints
	for (::Ipopt::Index i=0; i < boost::numeric_cast< ::Ipopt::Index>(m-m_pop->problem().get_ic_dimension());++i)
	{
		g_l[i] = g_u[i] = 0.0;
	}

	//Bounds on inequality constraints
	for (::Ipopt::Index i = ( m - m_pop->problem().get_ic_dimension()); i < m;++i)
	{
		g_l[i] = - 1e20;
		g_u[i] = 0.0;
	}

	return true;
}


bool ipopt_problem::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
					bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
					Ipopt::Index m, bool init_lambda,
					Ipopt::Number* lambda)
{
	// Here, we assume we only have starting values for x, if you code
	// your own NLP, you can provide starting values for the others if
	// you wish.
	pagmo_assert(init_x == true);
	pagmo_assert(init_z == false);
	pagmo_assert(init_lambda == false);
	pagmo_assert(n == boost::numeric_cast<Index>(m_pop->problem().get_dimension()));
	pagmo_assert(m == boost::numeric_cast<Index>(m_pop->problem().get_c_dimension()));

	// Shut off compiler warnings about unused variables. Cool, eh?
	(void) z_L;
	(void) z_U;
	(void) lambda;
	(void) n;
	(void) init_x;
	(void) init_z;
	(void) m;
	(void) init_lambda;

	// we initialize x as the best individual of the population
	pagmo::population::size_type bestidx = m_pop->get_best_idx();
	std::copy(m_pop->get_individual(bestidx).cur_x.begin(),m_pop->get_individual(bestidx).cur_x.end(),x);

	return true;
}


bool ipopt_problem::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value)
{
	(void) new_x;
	// return the value of the objective function
	std::copy(x,x+n,dv.begin());
	m_pop->problem().objfun(fit,dv);
	obj_value = fit[0];
	return true;
}

bool ipopt_problem::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f)
{
	(void) new_x;
	double central_diff;
	const double h0=1e-8;
	double h;
	std::copy(x,x+n,dv.begin());
	for (pagmo::decision_vector::size_type i=0; i<dv.size();++i)
	{
		grad_f[i] = 0;
	}

	double mem;
	for (size_t i =0;i<affects_obj.size();++i)
	{
		h = h0 * std::max(1.,fabs(dv[affects_obj[i]]));
		mem = dv[affects_obj[i]];
		dv[affects_obj[i]] += h;
		m_pop->problem().objfun(fit,dv);
		central_diff = fit[0];
		dv[affects_obj[i]] -= 2*h;
		m_pop->problem().objfun(fit,dv);
		central_diff = (central_diff-fit[0]) / 2 / h;
		grad_f[affects_obj[i]] = central_diff;
		dv[affects_obj[i]] = mem;
	}

	return true;
}

bool ipopt_problem::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g)
{
	(void) new_x;
	(void) m;
	// return the value of the constraints: g(x)
	std::copy(x,x+n,dv.begin());
	m_pop->problem().compute_constraints(con,dv);
	std::copy(con.begin(),con.end(),g);
	return true;
}

bool ipopt_problem::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
				Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
				Ipopt::Number* values)
{
	(void) new_x;
	(void) m;
	pagmo_assert(n == boost::numeric_cast<Index>(m_pop->problem().get_dimension()));
	pagmo_assert(m == boost::numeric_cast<Index>(m_pop->problem().get_c_dimension()));
	if (values == NULL) {
		// return the structure of the jacobian of the constraints
		for (Ipopt::Index i=0;i<nele_jac;++i)
		{
			iRow[i] = iJfun[i];
			jCol[i] = jJvar[i];
		}
	}
	else {
		double central_diff;
		const double h0 = 1e-8;
		double h;
		double mem;
		std::copy(x,x+n,dv.begin());
		for (Ipopt::Index i=0;i<nele_jac;++i)
		{
			h = h0 * std::max(1.,fabs(dv[jJvar[i]]));
			mem = dv[jJvar[i]];
			dv[jJvar[i]] += h;
			m_pop->problem().compute_constraints(con,dv);
			central_diff = con[iJfun[i]];
			dv[jJvar[i]] -= 2 * h;
			m_pop->problem().compute_constraints(con,dv);
			central_diff = (central_diff - con[iJfun[i]])/2/h;
			values[i] = central_diff;
			dv[jJvar[i]] = mem;
		}
	}

	return true;
}

void ipopt_problem::finalize_solution(Ipopt::SolverReturn status,
				       Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U,
				       Ipopt::Index m, const Ipopt::Number* g, const Ipopt::Number* lambda,
				       Ipopt::Number obj_value,
				       const Ipopt::IpoptData* ip_data,
				       Ipopt::IpoptCalculatedQuantities* ip_cq)
{
	(void) status;
	(void) z_L;
	(void) z_U;
	(void) m;
	(void) g;
	(void) lambda;
	(void) obj_value;
	(void) ip_data;
	(void) ip_cq;

	for (::Ipopt::Index i=0;i<n;i++) dv[i] = x[i];
	int bestidx = m_pop->get_best_idx();
	pagmo::decision_vector newx = dv;
	std::transform(dv.begin(), dv.end(), m_pop->get_individual(bestidx).cur_x.begin(), dv.begin(),std::minus<double>());
	for (int i=0;i<n;i++)
	{
		newx[i] = std::min(std::max(m_pop->problem().get_lb()[i],newx[i]),m_pop->problem().get_ub()[i]);
	}
	m_pop->set_x(bestidx,newx);
	m_pop->set_v(bestidx,dv);
}
