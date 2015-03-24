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

#include <boost/math/constants/constants.hpp>
#include <boost/lexical_cast.hpp>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "cec2009.h"

#define MYSIGN(x) ((x)>0?1.0:-1.0)

namespace pagmo { namespace problem {

static const double PI = boost::math::constants::pi<double>();

/// Constructor
/**
 * Will construct one of the 20 multi-objective optimization problems from
 * the CEC2009 competition. There are two sets of problems, namely the set
 * with unconstrained problems (UF) and the set with constrained problems (CF).
 *
 * @param[in] fun_id The problem id. One of [1,2,...10]
 * @param[in] prob_dim problem dimension. Default is 30, which is the setting used by the competition. But all the problems are scalable in terms of decision variable's dimension.
 * @param[in] is_constrained Specify whether the problem is constrained. False will yield the UF problems, and vice versa.
 *
 * @see http://www3.ntu.edu.sg/home/EPNSugan/index_files/CEC09-MOEA/CEC09-MOEA.htm
 *
 */
cec2009::cec2009(unsigned int fun_id, size_type prob_dim, bool is_constrained):
	base(prob_dim, 0, cec2009_fitness_dimension(fun_id), is_constrained?cec2009_ic_dimension(fun_id):0, is_constrained?cec2009_ic_dimension(fun_id):0, 0.0),
	m_problem_number(fun_id),
	m_is_constrained(is_constrained)
{
	configure_bounds();
}

/// Clone method
base_ptr cec2009::clone() const
{
	return base_ptr(new cec2009(*this));
}

std::string cec2009::get_name() const
{
	std::string retval("CEC2009 - ");
	if(!m_is_constrained){
		retval.append("UF");
	}
	else{
		retval.append("CF");
	}
	retval.append(boost::lexical_cast<std::string>(m_problem_number));
	return retval;
}

/// Returns the dimension of the fitness vector
fitness_vector::size_type cec2009::cec2009_fitness_dimension(int problem_id)
{
	// UF/CF 1-7 2-objective, UF/CF8-10 3-objective
	if(problem_id >= 1 && problem_id <= 7){
		return 2;
	}
	else if(problem_id >= 8 && problem_id <= 10){
		return 3;
	}
	else{
		pagmo_throw(value_error, "Error: CEC2009 unexpected problem id when determining fitness dimension.");
	}
}

/// Returns the number of inequality constraints for a particular CF problem
constraint_vector::size_type cec2009::cec2009_ic_dimension(int problem_id)
{
	// Note: Only come here if you are CF!
	// Only CF6 & CF7 have two inequality constraints, the rest has 1.
	// All the problems have no equality constraints.
	if(problem_id == 6 || problem_id == 7){
		return 2;
	}
	else if(problem_id >= 1 && problem_id <= 10){
		return 1;
	}
	else{
		pagmo_throw(value_error, "Error: CEC2009 unexpected problem id when determining inequality constraints\' dimension.");
	}
}

/// Configure the box bounds of the variables
void cec2009::configure_bounds(){

	size_type nx = get_dimension();
	std::vector<double> lb(nx, 0), ub(nx, 0);

	if(!m_is_constrained){
		if(m_problem_number == 1 || m_problem_number == 2 || m_problem_number == 5 || m_problem_number == 6 || m_problem_number == 7){
			// [0,1] x [-1,1]^{n-1}
			lb[0] = 0.0;
			ub[0] = 1.0;
			for(unsigned int i = 1; i < nx; i++){
				lb[i] = -1.0;
				ub[i] = 1.0;
			}
		}
		else if(m_problem_number == 3){
			// [0,1]^{n}
			for(unsigned int i = 0; i < nx; i++){
				lb[i] = 0.0;
				ub[i] = 1.0;
			}
		}
		else if(m_problem_number == 4){
			// [0,1] x [-2,2]^{n-1}
			lb[0] = 0.0;
			ub[0] = 1.0;
			for(unsigned int i = 1; i < nx; i++){
				lb[i] = -2.0;
				ub[i] = 2.0;
			}
		}
		else if(m_problem_number == 8 || m_problem_number == 9 || m_problem_number == 10){
			// [0,1]^{2} x [-2,2]^{n-2}
			lb[0] = 0.0;
			ub[0] = 1.0;
			lb[1] = 0.0;
			ub[1] = 1.0;
			for(unsigned int i = 2; i < nx; i++){
				lb[i] = -2.0;
				ub[i] = 2.0;
			}
		}
		else{
			pagmo_throw(value_error, "Error: CEC2009 unexpected problem id when setting bounds.");
		}
	}
	else{ // For CF
		if(m_problem_number == 2){
			// [0,1] x [-1,1]^{n-1}
			lb[0] = 0.0;
			ub[0] = 1.0;
			for(unsigned int i = 1; i < nx; i++){
				lb[i] = -1.0;
				ub[i] = 1.0;
			}
		}
		else if(m_problem_number == 1){
			// [0,1]^{n}
			for(unsigned int i = 0; i < nx; i++){
				lb[i] = 0.0;
				ub[i] = 1.0;
			}
		}
		else if(m_problem_number == 3 || m_problem_number == 4 || m_problem_number == 5 || m_problem_number == 6 || m_problem_number == 7){
			// [0,1] x [-2,2]^{n-1}
			lb[0] = 0.0;
			ub[0] = 1.0;
			for(unsigned int i = 1; i < nx; i++){
				lb[i] = -2.0;
				ub[i] = 2.0;
			}
		}
		else if(m_problem_number == 8){
			// [0,1]^{2} x [-4,4]^{n-2}
			lb[0] = 0.0;
			ub[0] = 1.0;
			lb[1] = 0.0;
			ub[1] = 1.0;
			for(unsigned int i = 2; i < nx; i++){
				lb[i] = -4.0;
				ub[i] = 4.0;
			}
		}
		else if(m_problem_number == 9 || m_problem_number == 10){
			// [0,1]^{2} x [-2,2]^{n-2}
			lb[0] = 0.0;
			ub[0] = 1.0;
			lb[1] = 0.0;
			ub[1] = 1.0;
			for(unsigned int i = 2; i < nx; i++){
				lb[i] = -2.0;
				ub[i] = 2.0;
			}
		}
		else{
			pagmo_throw(value_error, "Error: CEC2009 unexpected problem id when setting bounds.");
		}
	}
	set_bounds(lb, ub);
}

/// Implementation of the objective function.
void cec2009::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	size_type nx = get_dimension();
	if(!m_is_constrained){
		switch(m_problem_number)
		{
		case 1:
			UF1(&x[0], &f[0], nx);
			break;
		case 2:
			UF2(&x[0], &f[0], nx);
			break;
		case 3:
			UF3(&x[0], &f[0], nx);
			break;
		case 4:
			UF4(&x[0], &f[0], nx);
			break;
		case 5:
			UF5(&x[0], &f[0], nx);
			break;
		case 6:
			UF6(&x[0], &f[0], nx);
			break;
		case 7:
			UF7(&x[0], &f[0], nx);
			break;
		case 8:
			UF8(&x[0], &f[0], nx);
			break;
		case 9:
			UF9(&x[0], &f[0], nx);
			break;
		case 10:
			UF10(&x[0], &f[0], nx);
			break;
		default:
			pagmo_throw(value_error, "Error: CEC2009 unexpected problem id when computing fitness values.");
		}
	}
	else{
		switch(m_problem_number)
		{
		// NULL means do not compute the constraint values.
		case 1:
			CF1(&x[0], &f[0], NULL, nx);
			break;
		case 2:
			CF2(&x[0], &f[0], NULL, nx);
			break;
		case 3:
			CF3(&x[0], &f[0], NULL, nx);
			break;
		case 4:
			CF4(&x[0], &f[0], NULL, nx);
			break;
		case 5:
			CF5(&x[0], &f[0], NULL, nx);
			break;
		case 6:
			CF6(&x[0], &f[0], NULL, nx);
			break;
		case 7:
			CF7(&x[0], &f[0], NULL, nx);
			break;
		case 8:
			CF8(&x[0], &f[0], NULL, nx);
			break;
		case 9:
			CF9(&x[0], &f[0], NULL, nx);
			break;
		case 10:
			CF10(&x[0], &f[0], NULL, nx);
			break;
		default:
			pagmo_throw(value_error, "Error: CEC2009 unexpected problem id when computing fitness values.");
		}
	}
}

/// Implementation of the constraint computation
void cec2009::compute_constraints_impl(constraint_vector &c, const decision_vector & x) const
{
	size_type nx = get_dimension();
	if(!m_is_constrained){
		// TODO: Now just ignore; need to force c to zeros or throw?
		return;
	}
	// Allocation of f is done here as constraints computation requires the values of the objectives...
	// But here the process is forced to be de-coupled, hence recomputation of f is required.
	decision_vector f(nx);
	switch(m_problem_number)
	{
	case 1:
		CF1(&x[0], &f[0], &c[0], nx);
		break;
	case 2:
		CF2(&x[0], &f[0], &c[0], nx);
		break;
	case 3:
		CF3(&x[0], &f[0], &c[0], nx);
		break;
	case 4:
		CF4(&x[0], &f[0], &c[0], nx);
		break;
	case 5:
		CF5(&x[0], &f[0], &c[0], nx);
		break;
	case 6:
		CF6(&x[0], &f[0], &c[0], nx);
		break;
	case 7:
		CF7(&x[0], &f[0], &c[0], nx);
		break;
	case 8:
		CF8(&x[0], &f[0], &c[0], nx);
		break;
	case 9:
		CF9(&x[0], &f[0], &c[0], nx);
		break;
	case 10:
		CF10(&x[0], &f[0], &c[0], nx);
		break;
	default:
		pagmo_throw(value_error, "Error: CEC2009 unexpected problem id when computing constraint values.");
	}
}

// Below is slightly adapted from the problem code provided by the CEC2009
// competition organizing committee, e.g. changing the format of the
// constraint to g(x) <= 0.

/****************************************************************************/
// unconstraint test instances
/****************************************************************************/
void cec2009::UF1(const double *x, double *f, const unsigned int nx) const
{
	unsigned int j, count1, count2;
	double sum1, sum2, yj;
	
	sum1   = sum2   = 0.0;
	count1 = count2 = 0;
	for(j = 2; j <= nx; j++)
	{
		yj = x[j-1] - sin(6.0*PI*x[0] + j*PI/nx);
		yj = yj * yj;
		if(j % 2 == 0)
		{
			sum2 += yj;
			count2++;
		}
		else
		{
			sum1 += yj;
			count1++;
		}
	}
	f[0] = x[0]				+ 2.0 * sum1 / (double)count1;
	f[1] = 1.0 - sqrt(x[0]) + 2.0 * sum2 / (double)count2;
}

void cec2009::UF2(const double *x, double *f, const unsigned int nx) const
{
	unsigned int j, count1, count2;
	double sum1, sum2, yj;
	
	sum1   = sum2   = 0.0;
	count1 = count2 = 0;
	for(j = 2; j <= nx; j++)
	{
		if(j % 2 == 0) {
			yj = x[j-1]-0.3*x[0]*(x[0]*cos(24.0*PI*x[0]+4.0*j*PI/nx)+2.0)*sin(6.0*PI*x[0]+j*PI/nx);
			sum2 += yj*yj;
			count2++;
		}
		else
		{
			yj = x[j-1]-0.3*x[0]*(x[0]*cos(24.0*PI*x[0]+4.0*j*PI/nx)+2.0)*cos(6.0*PI*x[0]+j*PI/nx);
			sum1 += yj*yj;
			count1++;
		}
	}
	f[0] = x[0]				+ 2.0 * sum1 / (double)count1;
	f[1] = 1.0 - sqrt(x[0]) + 2.0 * sum2 / (double)count2;
}

void cec2009::UF3(const double *x, double *f, const unsigned int nx) const
{
	unsigned int j, count1, count2;
	double sum1, sum2, prod1, prod2, yj, pj;
	
	sum1   = sum2   = 0.0;
	count1 = count2 = 0;
	prod1  = prod2  = 1.0;
	for(j = 2; j <= nx; j++)
	{
		yj = x[j-1]-pow(x[0],0.5*(1.0+3.0*(j-2.0)/(nx-2.0)));
		pj = cos(20.0*yj*PI/sqrt(j+0.0));
		if (j % 2 == 0)
		{
			sum2  += yj*yj;
			prod2 *= pj;
			count2++;
		}
		else
		{
			sum1  += yj*yj;
			prod1 *= pj;
			count1++;
		}
	}
	f[0] = x[0]				+ 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1; f[1] = 1.0 - sqrt(x[0]) + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
}

void cec2009::UF4(const double *x, double *f, const unsigned int nx) const
{
	unsigned int j, count1, count2;
	double sum1, sum2, yj, hj;
	
	sum1   = sum2   = 0.0;
	count1 = count2 = 0;
	for(j = 2; j <= nx; j++)
	{
		yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
		hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
		if (j % 2 == 0)
		{
			sum2  += hj;
			count2++;
		}
		else
		{
			sum1  += hj;
			count1++;
		}
	}
	f[0] = x[0]				+ 2.0*sum1 / (double)count1;
	f[1] = 1.0 - x[0]*x[0]	+ 2.0*sum2 / (double)count2;
}

void cec2009::UF5(const double *x, double *f, const unsigned int nx) const
{
	unsigned int j, count1, count2;
	double sum1, sum2, yj, hj, N, E;
	
	sum1   = sum2   = 0.0;
	count1 = count2 = 0;
	N = 10.0; E = 0.1;
	for(j = 2; j <= nx; j++)
	{
		yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
		hj = 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
		if (j % 2 == 0)
		{
			sum2  += hj;
			count2++;
		}
		else
		{
			sum1  += hj;
			count1++;
		}
	}
	hj = (0.5/N + E)*fabs(sin(2.0*N*PI*x[0]));
	f[0] = x[0]	      + hj + 2.0*sum1 / (double)count1;
	f[1] = 1.0 - x[0] + hj + 2.0*sum2 / (double)count2;
}

void cec2009::UF6(const double *x, double *f, const unsigned int nx) const
{
	unsigned int j, count1, count2;
	double sum1, sum2, prod1, prod2, yj, hj, pj, N, E;
	N = 2.0; E = 0.1;

	sum1   = sum2   = 0.0;
	count1 = count2 = 0;
	prod1  = prod2  = 1.0;
	for(j = 2; j <= nx; j++)
	{
		yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
		pj = cos(20.0*yj*PI/sqrt(j+0.0));
		if (j % 2 == 0)
		{
			sum2  += yj*yj;
			prod2 *= pj;
			count2++;
		}
		else
		{
			sum1  += yj*yj;
			prod1 *= pj;
			count1++;
		}
	}

	hj = 2.0*(0.5/N + E)*sin(2.0*N*PI*x[0]);
	if(hj<0.0) hj = 0.0;
	f[0] = x[0]	      + hj + 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
	f[1] = 1.0 - x[0] + hj + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
}

void cec2009::UF7(const double *x, double *f, const unsigned int nx) const
{
	unsigned int j, count1, count2;
	double sum1, sum2, yj;
	
	sum1   = sum2   = 0.0;
	count1 = count2 = 0;
	for(j = 2; j <= nx; j++)
	{
		yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
		if (j % 2 == 0)
		{
			sum2  += yj*yj;
			count2++;
		}
		else
		{
			sum1  += yj*yj;
			count1++;
		}
	}
	yj = pow(x[0],0.2);
	f[0] = yj	    + 2.0*sum1 / (double)count1;
	f[1] = 1.0 - yj + 2.0*sum2 / (double)count2;
}

void cec2009::UF8(const double *x, double *f, const unsigned int nx) const
{
	unsigned int j, count1, count2, count3;
	double sum1, sum2, sum3, yj;
	
	sum1   = sum2   = sum3   = 0.0;
	count1 = count2 = count3 = 0;
	for(j = 3; j <= nx; j++)
	{
		yj = x[j-1] - 2.0*x[1]*sin(2.0*PI*x[0]+j*PI/nx);
		if(j % 3 == 1)
		{
			sum1  += yj*yj;
			count1++;
		}
		else if(j % 3 == 2)
		{
			sum2  += yj*yj;
			count2++;
		}
		else
		{
			sum3  += yj*yj;
			count3++;
		}
	}
	f[0] = cos(0.5*PI*x[0])*cos(0.5*PI*x[1]) + 2.0*sum1 / (double)count1;
	f[1] = cos(0.5*PI*x[0])*sin(0.5*PI*x[1]) + 2.0*sum2 / (double)count2;
	f[2] = sin(0.5*PI*x[0])                  + 2.0*sum3 / (double)count3;
}

void cec2009::UF9(const double *x, double *f, const unsigned int nx) const
{
	unsigned int j, count1, count2, count3;
	double sum1, sum2, sum3, yj, E;
	
	E = 0.1;
	sum1   = sum2   = sum3   = 0.0;
	count1 = count2 = count3 = 0;
	for(j = 3; j <= nx; j++)
	{
		yj = x[j-1] - 2.0*x[1]*sin(2.0*PI*x[0]+j*PI/nx);
		if(j % 3 == 1)
		{
			sum1  += yj*yj;
			count1++;
		}
		else if(j % 3 == 2)
		{
			sum2  += yj*yj;
			count2++;
		}
		else
		{
			sum3  += yj*yj;
			count3++;
		}
	}
	yj = (1.0+E)*(1.0-4.0*(2.0*x[0]-1.0)*(2.0*x[0]-1.0));
	if(yj<0.0) yj = 0.0;
	f[0] = 0.5*(yj + 2*x[0])*x[1]		+ 2.0*sum1 / (double)count1;
	f[1] = 0.5*(yj - 2*x[0] + 2.0)*x[1] + 2.0*sum2 / (double)count2;
	f[2] = 1.0 - x[1]                   + 2.0*sum3 / (double)count3;
}

void cec2009::UF10(const double *x, double *f, const unsigned int nx) const
{
	unsigned int j, count1, count2, count3;
	double sum1, sum2, sum3, yj, hj;
	
	sum1   = sum2   = sum3   = 0.0;
	count1 = count2 = count3 = 0;
	for(j = 3; j <= nx; j++)
	{
		yj = x[j-1] - 2.0*x[1]*sin(2.0*PI*x[0]+j*PI/nx);
		hj = 4.0*yj*yj - cos(8.0*PI*yj) + 1.0;
		if(j % 3 == 1)
		{
			sum1  += hj;
			count1++;
		}
		else if(j % 3 == 2)
		{
			sum2  += hj;
			count2++;
		}
		else
		{
			sum3  += hj;
			count3++;
		}
	}
	f[0] = cos(0.5*PI*x[0])*cos(0.5*PI*x[1]) + 2.0*sum1 / (double)count1;
	f[1] = cos(0.5*PI*x[0])*sin(0.5*PI*x[1]) + 2.0*sum2 / (double)count2;
	f[2] = sin(0.5*PI*x[0])                  + 2.0*sum3 / (double)count3;
}

/****************************************************************************/
// constraint test instances
/****************************************************************************/
void cec2009::CF1(const double *x, double *f, double *c, const unsigned int nx) const
{
	unsigned int j, count1, count2;
	double sum1, sum2, yj, N, a;
	N = 10.0; a = 1.0;
	
	sum1   = sum2   = 0.0;
	count1 = count2 = 0;
	for(j = 2; j <= nx; j++)
	{
		yj = x[j-1]-pow(x[0],0.5*(1.0+3.0*(j-2.0)/(nx-2.0)));
		if (j % 2 == 1)
		{
			sum1  += yj*yj;
			count1++;
		}
		else
		{
			sum2  += yj*yj;
			count2++;
		}
	}
	f[0] = x[0]		  + 2.0*sum1 / (double)count1;
	f[1] = 1.0 - x[0] + 2.0*sum2 / (double)count2;

	if(c != NULL){
		c[0] = f[1] + f[0] - a*fabs(sin(N*PI*(f[0]-f[1]+1.0))) - 1.0; 
		c[0] = -c[0]; //convert to g(x) <= 0 form
	}
}

void cec2009::CF2(const double *x, double *f, double *c, const unsigned int nx) const
{
	unsigned int j, count1, count2;
	double sum1, sum2, yj, N, a, t;
	N = 2.0; a = 1.0;
	
	sum1   = sum2   = 0.0;
	count1 = count2 = 0;
	for(j = 2; j <= nx; j++)
	{
		if (j % 2 == 1)
		{
			yj = x[j-1] - sin(6.0*PI*x[0] + j*PI/nx);
			sum1  += yj*yj;
			count1++;
		}
		else
		{
			yj = x[j-1] - cos(6.0*PI*x[0] + j*PI/nx);
			sum2  += yj*yj;
			count2++;
		}
	}
	f[0] = x[0]		        + 2.0*sum1 / (double)count1;
	f[1] = 1.0 - sqrt(x[0]) + 2.0*sum2 / (double)count2;

	if(c != NULL){
		t	 = f[1] + sqrt(f[0]) - a*sin(N*PI*(sqrt(f[0])-f[1]+1.0)) - 1.0;
		c[0] = MYSIGN(t)*fabs(t)/(1+exp(4.0*fabs(t)));
		c[0] = -c[0]; //convert to g(x) <= 0 form
	}
}

void cec2009::CF3(const double *x, double *f, double *c, const unsigned int nx) const
{
	unsigned int j, count1, count2;
	double sum1, sum2, prod1, prod2, yj, pj, N, a;
	N = 2.0; a = 1.0;

	sum1   = sum2   = 0.0;
	count1 = count2 = 0;
	prod1  = prod2  = 1.0;
	for(j = 2; j <= nx; j++)
	{
		yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
		pj = cos(20.0*yj*PI/sqrt(j+0.0));
		if (j % 2 == 0)
		{
			sum2  += yj*yj;
			prod2 *= pj;
			count2++;
		}
		else
		{
			sum1  += yj*yj;
			prod1 *= pj;
			count1++;
		}
	}

	f[0] = x[0]	           + 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
	f[1] = 1.0 - x[0]*x[0] + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;

	if(c != NULL){
		c[0] = f[1] + f[0]*f[0] - a*sin(N*PI*(f[0]*f[0]-f[1]+1.0)) - 1.0;
		c[0] = -c[0]; //convert to g(x) <= 0 form
	}
}

void cec2009::CF4(const double *x, double *f, double *c, const unsigned int nx) const
{
	unsigned int j;
	double sum1, sum2, yj, t;

	sum1   = sum2   = 0.0;
	for(j = 2; j <= nx; j++)
	{
		yj = x[j-1] - sin(6.0*PI*x[0] + j*PI/nx);
		if (j % 2 == 1)
		{
			sum1  += yj*yj;
		}
		else
		{
			if (j==2)
				sum2 += yj < 1.5-0.75*sqrt(2.0) ? fabs(yj) : (0.125+(yj-1)*(yj-1));
			else
				sum2  += yj*yj;
		}
	}
	f[0] = x[0]		  + sum1;
	f[1] = 1.0 - x[0] + sum2;

	if(c != NULL){
		t	 = x[1] - sin(6.0*x[0]*PI+2.0*PI/nx) - 0.5*x[0] + 0.25;
		c[0] = MYSIGN(t)*fabs(t)/(1+exp(4.0*fabs(t)));
		c[0] = -c[0]; //convert to g(x) <= 0 form
	}
}

void cec2009::CF5(const double *x, double *f, double *c, const unsigned int nx) const
{
	unsigned int j;
	double sum1, sum2, yj;

	sum1   = sum2   = 0.0;
	for(j = 2; j <= nx; j++)
	{
		if (j % 2 == 1)
		{
			yj    = x[j-1] - 0.8*x[0]*cos(6.0*PI*x[0] + j*PI/nx);
			sum1 += 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
		}
		else
		{
			yj = x[j-1] - 0.8*x[0]*sin(6.0*PI*x[0] + j*PI/nx);
			if (j==2)
				sum2 += yj < 1.5-0.75*sqrt(2.0) ? fabs(yj) : (0.125+(yj-1)*(yj-1));
			else
				sum2 += 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
		}
	}
	f[0] = x[0]		  + sum1;
	f[1] = 1.0 - x[0] + sum2;

	if(c != NULL){
		c[0] = x[1] - 0.8*x[0]*sin(6.0*x[0]*PI+2.0*PI/nx) - 0.5*x[0] + 0.25;
		c[0] = -c[0]; //convert to g(x) <= 0 form
	}
}

void cec2009::CF6(const double *x, double *f, double *c, const unsigned int nx) const
{
	unsigned int j;
	double sum1, sum2, yj;
	
	sum1   = sum2   = 0.0;
	for(j = 2; j <= nx; j++)
	{
		if (j % 2 == 1)
		{
			yj     = x[j-1] - 0.8*x[0]*cos(6.0*PI*x[0] + j*PI/nx);
			sum1  += yj*yj;
		}
		else
		{
			yj     = x[j-1] - 0.8*x[0]*sin(6.0*PI*x[0] + j*PI/nx);
			sum2  += yj*yj;
		}
	}
	f[0] = x[0]		                 + sum1;
	f[1] = (1.0 - x[0])*(1.0 - x[0]) + sum2;

	if(c != NULL){
		c[0] = x[1]-0.8*x[0]*sin(6.0*x[0]*PI+2.0*PI/nx) - MYSIGN((x[0]-0.5)*(1.0-x[0]))*sqrt(fabs((x[0]-0.5)*(1.0-x[0])));
		c[1] = x[3]-0.8*x[0]*sin(6.0*x[0]*PI+4.0*PI/nx) - MYSIGN(0.25*sqrt(1-x[0])-0.5*(1.0-x[0]))*sqrt(fabs(0.25*sqrt(1-x[0])-0.5*(1.0-x[0])));
		//convert to g(x) <= 0 form
		c[0] = -c[0];
		c[1] = -c[1];
	}
}

void cec2009::CF7(const double *x, double *f, double *c, const unsigned int nx) const
{
	unsigned int j;
	double sum1, sum2, yj;

	sum1   = sum2   = 0.0;
	for(j = 2; j <= nx; j++)
	{
		if (j % 2 == 1)
		{
			yj     = x[j-1] - cos(6.0*PI*x[0] + j*PI/nx);
			sum1  += 2.0*yj*yj-cos(4.0*PI*yj)+1.0;
		}
		else
		{
			yj     = x[j-1] - sin(6.0*PI*x[0] + j*PI/nx);
			if (j==2 || j==4)
				sum2 += yj*yj;
			else
				sum2  += 2.0*yj*yj-cos(4.0*PI*yj)+1.0;
		}
	}
	f[0] = x[0]		                 + sum1;
	f[1] = (1.0 - x[0])*(1.0 - x[0]) + sum2;

	if(c != NULL){
		c[0] = x[1]-sin(6.0*x[0]*PI+2.0*PI/nx) - MYSIGN((x[0]-0.5)*(1.0-x[0]))*sqrt(fabs((x[0]-0.5)*(1.0-x[0])));
		c[1] = x[3]-sin(6.0*x[0]*PI+4.0*PI/nx) - MYSIGN(0.25*sqrt(1-x[0])-0.5*(1.0-x[0]))*sqrt(fabs(0.25*sqrt(1-x[0])-0.5*(1.0-x[0])));
		//convert to g(x) <= 0 form
		c[0] = -c[0];
		c[1] = -c[1];
	}
}

void cec2009::CF8(const double *x, double *f, double *c, const unsigned int nx) const
{
	unsigned int j, count1, count2, count3;
	double sum1, sum2, sum3, yj, N, a;
	N = 2.0; a = 4.0;

	sum1   = sum2   = sum3   = 0.0;
	count1 = count2 = count3 = 0;
	for(j = 3; j <= nx; j++)
	{
		yj = x[j-1] - 2.0*x[1]*sin(2.0*PI*x[0]+j*PI/nx);
		if(j % 3 == 1)
		{
			sum1  += yj*yj;
			count1++;
		}
		else if(j % 3 == 2)
		{
			sum2  += yj*yj;
			count2++;
		}
		else
		{
			sum3  += yj*yj;
			count3++;
		}
	}
	f[0] = cos(0.5*PI*x[0])*cos(0.5*PI*x[1]) + 2.0*sum1 / (double)count1;
	f[1] = cos(0.5*PI*x[0])*sin(0.5*PI*x[1]) + 2.0*sum2 / (double)count2;
	f[2] = sin(0.5*PI*x[0])                  + 2.0*sum3 / (double)count3;

	if(c != NULL){
		c[0] = (f[0]*f[0]+f[1]*f[1])/(1-f[2]*f[2]) - a*fabs(sin(N*PI*((f[0]*f[0]-f[1]*f[1])/(1-f[2]*f[2])+1.0))) - 1.0;
		c[0] = -c[0]; //convert to g(x) <= 0 form
	}
}

void cec2009::CF9(const double *x, double *f, double *c, const unsigned int nx) const
{
	unsigned int j, count1, count2, count3;
	double sum1, sum2, sum3, yj, N, a;
	N = 2.0; a = 3.0;

	sum1   = sum2   = sum3   = 0.0;
	count1 = count2 = count3 = 0;
	for(j = 3; j <= nx; j++)
	{
		yj = x[j-1] - 2.0*x[1]*sin(2.0*PI*x[0]+j*PI/nx);
		if(j % 3 == 1)
		{
			sum1  += yj*yj;
			count1++;
		}
		else if(j % 3 == 2)
		{
			sum2  += yj*yj;
			count2++;
		}
		else
		{
			sum3  += yj*yj;
			count3++;
		}
	}
	f[0] = cos(0.5*PI*x[0])*cos(0.5*PI*x[1]) + 2.0*sum1 / (double)count1;
	f[1] = cos(0.5*PI*x[0])*sin(0.5*PI*x[1]) + 2.0*sum2 / (double)count2;
	f[2] = sin(0.5*PI*x[0])                  + 2.0*sum3 / (double)count3;

	if(c != NULL){
		c[0] = (f[0]*f[0]+f[1]*f[1])/(1-f[2]*f[2]) - a*sin(N*PI*((f[0]*f[0]-f[1]*f[1])/(1-f[2]*f[2])+1.0)) - 1.0;
		c[0] = -c[0]; //convert to g(x) <= 0 form
	}
}

void cec2009::CF10(const double *x, double *f, double *c, const unsigned int nx) const
{
	unsigned int j, count1, count2, count3;
	double sum1, sum2, sum3, yj, hj, N, a;
	N = 2.0; a = 1.0;

	sum1   = sum2   = sum3   = 0.0;
	count1 = count2 = count3 = 0;
	for(j = 3; j <= nx; j++)
	{
		yj = x[j-1] - 2.0*x[1]*sin(2.0*PI*x[0]+j*PI/nx);
		hj = 4.0*yj*yj - cos(8.0*PI*yj) + 1.0;
		if(j % 3 == 1)
		{
			sum1  += hj;
			count1++;
		}
		else if(j % 3 == 2)
		{
			sum2  += hj;
			count2++;
		}
		else
		{
			sum3  += hj;
			count3++;
		}
	}
	f[0] = cos(0.5*PI*x[0])*cos(0.5*PI*x[1]) + 2.0*sum1 / (double)count1;
	f[1] = cos(0.5*PI*x[0])*sin(0.5*PI*x[1]) + 2.0*sum2 / (double)count2;
	f[2] = sin(0.5*PI*x[0])                  + 2.0*sum3 / (double)count3;
	
	if(c != NULL){
		c[0] = (f[0]*f[0]+f[1]*f[1])/(1-f[2]*f[2]) - a*sin(N*PI*((f[0]*f[0]-f[1]*f[1])/(1-f[2]*f[2])+1.0)) - 1.0;
		c[0] = -c[0]; //convert to g(x) <= 0 form
	}
}

}} //namespaces

#undef MYSIGN

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::cec2009)
