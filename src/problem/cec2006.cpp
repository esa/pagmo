/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
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

#include <cmath>

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "cec2006.h"

namespace pagmo { namespace problem {

const double PI = 4.* std::atan(1.);

// initialization of the problems, global constraints and inequality constraints dimension
const decision_vector::size_type cec2006::m_problems_dimension[] =
{13,20,10,5,4,2,10,2,7,8,2,3,5,10,3,5,6,9,15,24,7,22,9,2};
const constraint_vector::size_type cec2006::m_problems_c_dimension[] =
{9,2,1,6,5,2,8,2,4,6,1,1,3,3,2,38,4,13,5,20,6,20,6,2};
const constraint_vector::size_type cec2006::m_problems_ic_dimension[] =
{9,2,0,6,2,2,8,2,4,6,0,1,0,0,0,38,0,13,5,6,1,1,2,2};

/// Constructor
/**
 * Will construct one of the 24 CEC2006 problems
 *
 * @param[in] fun_id The problem id. One of [1,2,...,24]
 */
cec2006::cec2006(unsigned int fun_id):base(m_problems_dimension[fun_id],0,1,m_problems_c_dimension[fun_id],m_problems_ic_dimension[fun_id]),m_problem_number(fun_id)
{
    // set the bounds for the current problem
    switch(m_problem_number)
    {
    case 1:
    {
        const double lb[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
        const double ub[] = {1.,1.,1.,1.,1.,1.,1.,1.,1.,100.,100.,100.,1.};
        this->set_bounds(lb,ub);
        break;
    }
    case 2:
    {
        const double lb[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
        const double ub[] = {10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.};
        this->set_bounds(lb,ub);
        break;
    }
    case 3:
    {
        const double lb[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
        const double ub[] = {1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};
        this->set_bounds(lb,ub);
        break;
    }
    case 4:
    {
        const double lb[] = {78.,33.,27.,27.,27.};
        const double ub[] = {102.,45.,45.,45.,45.};
        this->set_bounds(lb,ub);
        break;
    }
    case 5:
    {
        const double lb[] = {0.,0.,-0.55,-0.55};
        const double ub[] = {1.,1.,1.,1.};
        this->set_bounds(lb,ub);
        break;
    }
    case 6:
    {
        const double lb[] = {13.,0.};
        const double ub[] = {100.,100.};
        this->set_bounds(lb,ub);
        break;
    }
    case 7:
    {
        const double lb[] = {-10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.};
        const double ub[] = {10.,10.,10.,10.,10.,10.,10.,10.,10.,10.};
        this->set_bounds(lb,ub);
        break;
    }
    case 8:
    {
        const double lb[] = {0.,0.};
        const double ub[] = {10.,10.};
        this->set_bounds(lb,ub);
        break;
    }
    case 9:
    {
        const double lb[] = {-10.,-10.,-10.,-10.,-10.,-10.,-10.};
        const double ub[] = {10.,10.,10.,10.,10.,10.,10.};
        this->set_bounds(lb,ub);
        break;
    }
    case 10:
    {
        const double lb[] = {100.,1000.,1000.,10.,10.,10.,10.,10.};
        const double ub[] = {10000.,10000.,10000.,1000.,1000.,1000.,1000.,1000.};
        this->set_bounds(lb,ub);
        break;
    }
    case 11:
    {
        const double lb[] = {-1.,-1.};
        const double ub[] = {1.,1.};
        this->set_bounds(lb,ub);
        break;
    }
    case 12:
    {
        const double lb[] = {0.,0.,0.};
        const double ub[] = {10.,10.,10.};
        this->set_bounds(lb,ub);
        break;
    }
    case 13:
    {
        const double lb[] = {-2.3,-2.3,-3.2,-3.2,-3.2};
        const double ub[] = {2.3,2.3,3.2,3.2,3.2};
        this->set_bounds(lb,ub);
        break;
    }
    case 14:
    {
        const double lb[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
        const double ub[] = {10.,10.,10.,10.,10.,10.,10.,10.,10.,10.};
        this->set_bounds(lb,ub);
        break;
    }
    case 15:
    {
        const double lb[] = {0.,0.,0.};
        const double ub[] = {10.,10.,10.};
        this->set_bounds(lb,ub);
        break;
    }
    case 16:
    {
        const double lb[] = {704.4148,68.6,0.,193.,25.};
        const double ub[] = {906.3855,288.88,134.75,287.0966,84.1988};
        this->set_bounds(lb,ub);
        break;
    }
    default:
        pagmo_throw(value_error, "Error: There are only 24 test functions in this test suite!");
        break;
    }
}

/// Clone method.
base_ptr cec2006::clone() const
{
    return base_ptr(new cec2006(*this));
}

/// Implementation of the objective function.
void cec2006::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
    size_type nx = get_dimension();
    switch(m_problem_number)
    {
    case 1:
        this->g01_objfun_impl(f,x);
        break;
    case 2:
        this->g02_objfun_impl(f,x);
        break;
    case 3:
        this->g03_objfun_impl(f,x);
        break;
    case 4:
        this->g04_objfun_impl(f,x);
        break;
    default:
        pagmo_throw(value_error, "Error: There are only 24 test functions in this test suite!");
        break;
    }
}

/// Implementation of the constraint function.
void cec2006::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
    switch(m_problem_number)
    {
    case 1:
        this->g01_compute_constraints_impl(c,x);
        break;
    case 2:
        break;
    default:
        pagmo_throw(value_error, "Error: There are only 24 test functions in this test suite!");
        break;
    }
}

std::string cec2006::get_name() const
{
    std::string retval("CEC2006 - g");
    retval.append(boost::lexical_cast<std::string>(m_problem_number));

    return retval;
}

// -------------------------------------------

/// Implementation of the objective function.
void cec2006::g01_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
    /* objective function */
    f[0] = 5.0 * (x[0] + x[1] + x[2] + x[3]) - 5.0 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);

    for (decision_vector::size_type j = 4; j < 13; j++)
        f[0] = f[0] - x[j];
}

/// Implementation of the constraint function.
void cec2006::g01_compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
    /* constraints g<=0 */
    c[0] = 2.0 * x[0] + 2.0 * x[1] + x[9] + x[10] - 10.;
    c[1] = 2.0 * x[0] + 2.0 * x[2] + x[9] + x[11] - 10.;
    c[2] = 2.0 * x[1] + 2.0 * x[2] + x[10] + x[11] - 10.;
    c[3] = -8.0 * x[0] + x[9];
    c[4] = -8.0 * x[1] + x[10];
    c[5] = -8.0 * x[2] + x[11];
    c[6] = -2.0 * x[3] - x[4] + x[9];
    c[7] = -2.0 * x[5] - x[6] + x[10];
    c[8] = -2.0 * x[7] - x[8] + x[11];
}

// -------------------------------------------

/// Implementation of the objective function.
void cec2006::g02_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
    double nx = this->get_dimension();

    /* objective function */
    double f1 = 0.;
    double f2 = 1.;
    double f3 = 0.;

    for (int j = 0; j < nx; j++)
    {
        f1 = f1 + pow (cos (x[j]), 4);
        f2 = f2 * cos (x[j]) * cos (x[j]);
        f3 = f3 + ((double) (j + 1)) * x[j] * x[j];
    }
    f[0] = - fabs ((f1 - 2 * f2) / sqrt (f3));
}

/// Implementation of the constraint function.
void cec2006::g02_compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
    /* constraints g<=0 */
    double nx = this->get_dimension();

    double g1 = 1.;
    double g2 = 0.;

    for (int j = 0; j < nx; j++)
    {
        g1 = g1 * x[j];
        g2 = g2 + x[j];
    }

    c[0] = 0.75 - g1;
    c[1] = g2 - 7.5 * ((double) nx);
}

// -------------------------------------------

/// Implementation of the objective function.
void cec2006::g03_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
    /* objective function */
    double nx = this->get_dimension();

    double f1 = 1.;
    double f3 = sqrt ((double) nx);

    for (int j = 0; j < nx; j++)
    {
        f1 = f3 * f1 * x[j];
    }

    f[0] = - f1;
}

/// Implementation of the constraint function.
void cec2006::g03_compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
    double nx = this->get_dimension();

    double f2 = 0.;

    for (int j = 0; j < nx; j++)
    {
        f2 = f2 + x[j] * x[j];
    }

    /* constraints h=0 */
    c[0] = f2 - 1.0;
}

// -------------------------------------------

/// Implementation of the objective function.
void cec2006::g04_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
    /* objective function */
    f[0] = 5.3578547 * x[2] * x[2] + 0.8356891 * x[0] * x[4] + 37.293239 * x[0] - 40792.141;
}

/// Implementation of the constraint function.
void cec2006::g04_compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
    /* constraints g<=0 */
    c[0] = 85.334407 + 0.0056858 * x[1] * x[4] + 0.0006262 * x[0] * x[3] - 0.0022053 * x[2] * x[4] - 92.;
    c[1] = -85.334407 - 0.0056858 * x[1] * x[4] - 0.0006262 * x[0] * x[3] + 0.0022053 * x[2] * x[4];
    c[2] = 80.51249 + 0.0071317 * x[1] * x[4] + 0.0029955 * x[0] * x[1] + 0.0021813 * x[2] * x[2] - 110.;
    c[3] = -80.51249 - 0.0071317 * x[1] * x[4] - 0.0029955 * x[0] * x[1] - 0.0021813 * x[2] * x[2] + 90.;
    c[4] = 9.300961 + 0.0047026 * x[2] * x[4] + 0.0012547 * x[0] * x[2] + 0.0019085 * x[2] * x[3] - 25.;
    c[5] = -9.300961 - 0.0047026 * x[2] * x[4] - 0.0012547 * x[0] * x[2] - 0.0019085 * x[2] * x[3] + 20.;
}

// -------------------------------------------

/// Implementation of the objective function.
void cec2006::g05_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
    /* objective function */
    f[0] = 3.0 * x[0] + 0.000001 * pow (x[0], 3) + 2.0 * x[1] + (0.000002 / 3.0) * pow (x[1], 3);
}

/// Implementation of the constraint function.
void cec2006::g05_compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
    /* constraints h=0 */
    c[0] = 1000.0 * sin (-x[2] - 0.25) + 1000.0 * sin (-x[3] - 0.25) + 894.8 - x[0];
    c[1] = 1000.0 * sin (x[2] - 0.25) + 1000.0 * sin (x[2] - x[3] - 0.25) + 894.8 - x[1];
    c[2] = 1000.0 * sin (x[3] - 0.25) + 1000.0 * sin (x[3] - x[2] - 0.25) + 1294.8;

    /* constraints g<=0 */
    c[3] = -x[3] + x[2] - 0.55;
    c[4] = -x[2] + x[3] - 0.55;
}

// -------------------------------------------

/// Implementation of the objective function.
void cec2006::g06_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
    /* objective function */
    f[0] = pow ((x[0] - 10.), 3) + pow ((x[1] - 20.), 3);
}

/// Implementation of the constraint function.
void cec2006::g06_compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
    /* constraints g<=0 */
    c[0] = 100. - (x[0] - 5.) * (x[0] - 5.) - (x[1] - 5.) * (x[1] - 5.);
    c[1] = (x[0] - 6.) * (x[0] - 6.) + (x[1] - 5.) * (x[1] - 5.) - 82.81;
}

// -------------------------------------------

/// Implementation of the objective function.
void cec2006::g07_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
    /* objective function */
    f[0] =
            x[0] * x[0] + x[1] * x[1] + x[0] * x[1] - 14.0 * x[0] - 16.0 * x[1] + (x[2] - 10.0) * (x[2] - 10.0) + 4.0 * (x[3] -
            5.0) *
            (x[3] - 5.0) + (x[4] - 3.0) * (x[4] - 3.0) + 2.0 * (x[5] - 1.0) * (x[5] - 1.0) + 5.0 * x[6] * x[6] + 7.0 * (x[7] -
            11) *
            (x[7] - 11) + 2.0 * (x[8] - 10.0) * (x[8] - 10.0) + (x[9] - 7.0) * (x[9] - 7.0) + 45.;
}

/// Implementation of the constraint function.
void cec2006::g07_compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
    /* constraints g<=0 */
    c[0] = -105.0 + 4.0 * x[0] + 5.0 * x[1] - 3.0 * x[6] + 9.0 * x[7];
    c[1] = 10.0 * x[0] - 8.0 * x[1] - 17.0 * x[6] + 2.0 * x[7];
    c[2] = -8.0 * x[0] + 2.0 * x[1] + 5.0 * x[8] - 2.0 * x[9] - 12.0;
    c[3] = 3.0 * (x[0] - 2.0) * (x[0] - 2.0) + 4.0 * (x[1] - 3.0) * (x[1] - 3.0) + 2.0 * x[2] * x[2] - 7.0 * x[3] - 120.0;
    c[4] = 5.0 * x[0] * x[0] + 8.0 * x[1] + (x[2] - 6.0) * (x[2] - 6.0) - 2.0 * x[3] - 40.0;
    c[5] = x[0] * x[0] + 2.0 * (x[1] - 2.0) * (x[1] - 2.0) - 2.0 * x[0] * x[1] + 14.0 * x[4] - 6.0 * x[5];
    c[6] = 0.5 * (x[0] - 8.0) * (x[0] - 8.0) + 2.0 * (x[1] - 4.0) * (x[1] - 4.0) + 3.0 * x[4] * x[4] - x[5] - 30.0;
    c[7] = -3.0 * x[0] + 6.0 * x[1] + 12.0 * (x[8] - 8.0) * (x[8] - 8.0) - 7.0 * x[9];
}

// -------------------------------------------

/// Implementation of the objective function.
void cec2006::g08_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
    /* objective function */
    f[0] = - std::pow(std::sin (2 * PI * x[0]), 3) * std::sin (2 * PI * x[1]) / ( std::pow (x[0], 3) * (x[0] + x[1]) );
}

/// Implementation of the constraint function.
void cec2006::g08_compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
    /* constraints g<=0 */
    c[0] = x[0] * x[0] - x[1] + 1.0;
    c[1] = 1.0 - x[0] + (x[1] - 4.0) * (x[1] - 4.0);
}

//// -------------------------------------------

///// Default constructor
//cec2006_g09::cec2006_g09():base(7,0,1,4,4)
//{
//}

///// Clone method.
//base_ptr cec2006_g09::clone() const
//{
//    return base_ptr(new cec2006_g09(*this));
//}

///// Implementation of the objective function.
//void cec2006_g09::objfun_impl(fitness_vector &f, const decision_vector &x) const
//{
//    /* objective function */
//    f[0] = (x[0] - 10.0) * (x[0] - 10.0) + 5.0 * (x[1] - 12.0) * (x[1] - 12.0) + pow (x[2], 4) +
//            3.0 * (x[3] - 11.0) * (x[3] - 11.0) + 10.0 * pow (x[4], 6) + 7.0 * x[5] * x[5] +
//            pow (x[6], 4) - 4.0 * x[5] * x[6] - 10.0 * x[5] - 8.0 * x[6];
//}

///// Implementation of the constraint function.
//void cec2006_g09::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
//{
//    /* constraints g<=0 */
//    c[0] = -127.0 + 2 * x[0] * x[0] + 3.0 * pow (x[1], 4) + x[2] + 4.0 * x[3] * x[3] + 5.0 * x[4];
//    c[1] = -282.0 + 7.0 * x[0] + 3.0 * x[1] + 10.0 * x[2] * x[2] + x[3] - x[4];
//    c[2] = -196.0 + 23.0 * x[0] + x[1] * x[1] + 6.0 * x[5] * x[5] - 8.0 * x[6];
//    c[3] = 4.0 * x[0] * x[0] + x[1] * x[1] - 3.0 * x[0] * x[1] + 2.0 * x[2] * x[2] + 5.0 * x[5] - 11.0 * x[6];
//}

//std::string cec2006_g09::get_name() const
//{
//    return "CEC2006 - g09";
//}

//// -------------------------------------------

///// Default constructor
//cec2006_g10::cec2006_g10():base(8,0,1,6,6)
//{
//}

///// Clone method.
//base_ptr cec2006_g10::clone() const
//{
//    return base_ptr(new cec2006_g10(*this));
//}

///// Implementation of the objective function.
//void cec2006_g10::objfun_impl(fitness_vector &f, const decision_vector &x) const
//{
//    /* objective function */
//    f[0] = x[0] + x[1] + x[2];
//}

///// Implementation of the constraint function.
//void cec2006_g10::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
//{
//    /* constraints g<=0 */
//    c[0] = -1.0 + 0.0025 * (x[3] + x[5]);
//    c[1] = -1.0 + 0.0025 * (x[4] + x[6] - x[3]);
//    c[2] = -1.0 + 0.01 * (x[7] - x[4]);
//    c[3] = -x[0] * x[5] + 833.33252 * x[3] + 100.0 * x[0] - 83333.333;
//    c[4] = -x[1] * x[6] + 1250.0 * x[4] + x[1] * x[3] - 1250.0 * x[3];
//    c[5] = -x[2] * x[7] + 1250000.0 + x[2] * x[4] - 2500.0 * x[4];
//}

//std::string cec2006_g10::get_name() const
//{
//    return "CEC2006 - g10";
//}

//// -------------------------------------------

///// Default constructor
//cec2006_g11::cec2006_g11():base(2,0,1,1,0)
//{
//}

///// Clone method.
//base_ptr cec2006_g11::clone() const
//{
//    return base_ptr(new cec2006_g11(*this));
//}

///// Implementation of the objective function.
//void cec2006_g11::objfun_impl(fitness_vector &f, const decision_vector &x) const
//{
//    /* objective function */
//    f[0] = x[0] * x[0] + (x[1] - 1.0) * (x[1] - 1.0);
//}

///// Implementation of the constraint function.
//void cec2006_g11::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
//{
//    /* constraints h=0 */
//    c[0] = x[1] - x[0] * x[0];
//}

//std::string cec2006_g11::get_name() const
//{
//    return "CEC2006 - g11";
//}

//// -------------------------------------------

///// Default constructor
//cec2006_g12::cec2006_g12():base(3,0,1,1,1)
//{
//}

///// Clone method.
//base_ptr cec2006_g12::clone() const
//{
//    return base_ptr(new cec2006_g12(*this));
//}

///// Implementation of the objective function.
//void cec2006_g12::objfun_impl(fitness_vector &f, const decision_vector &x) const
//{
//    /* objective function */
//    f[0] = - (100. - (x[0] - 5.) * (x[0] - 5.) - (x[1] - 5.) * (x[1] - 5.) - (x[2] - 5.) * (x[2] - 5.)) / 100.;
//}

///// Implementation of the constraint function.
//void cec2006_g12::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
//{
//    double gt;

//    /* constraints g<=0 */
//    g[0] = (x[0] - 1.) * (x[0] - 1.) + (x[1] - 1.) * (x[1] - 1.) + (x[2] - 1.) * (x[2] - 1.) - 0.0625;
//    for (int i = 1; i <= 9; i++)
//    {
//        for (int j = 1; j <= 9; j++)
//        {
//            for (int k = 1; k <= 9; k++)
//            {
//                gt = (x[0] - i) * (x[0] - i) + (x[1] - j) * (x[1] - j) + (x[2] - k) * (x[2] - k) - 0.0625;
//                if (gt < g[0])
//                    g[0] = gt;
//            }
//        }
//    }
//}

//std::string cec2006_g12::get_name() const
//{
//    return "CEC2006 - g12";
//}

//// -------------------------------------------

///// Default constructor
//cec2006_g13::cec2006_g13():base(5,0,1,3,0)
//{
//}

///// Clone method.
//base_ptr cec2006_g13::clone() const
//{
//    return base_ptr(new cec2006_g13(*this));
//}

///// Implementation of the objective function.
//void cec2006_g13::objfun_impl(fitness_vector &f, const decision_vector &x) const
//{
//    /* objective function */
//    f[0] = exp(x[0] * x[1] * x[2] * x[3] * x[4]);
//}

///// Implementation of the constraint function.
//void cec2006_g13::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
//{
//    /* constraints h(x) = 0 */
//    c[0] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + x[4] * x[4] - 10.0;
//    c[1] = x[1] * x[2] - 5.0 * x[3] * x[4];
//    c[2] = pow(x[0], 3) + pow(x[1], 3) + 1.0;
//}

//std::string cec2006_g13::get_name() const
//{
//    return "CEC2006 - g13";
//}

//// -------------------------------------------

///// Default constructor
//cec2006_g14::cec2006_g14():base(10,0,1,3,0)
//{
//}

///// Clone method.
//base_ptr cec2006_g14::clone() const
//{
//    return base_ptr(new cec2006_g14(*this));
//}

///// Implementation of the objective function.
//void cec2006_g14::objfun_impl(fitness_vector &f, const decision_vector &x) const
//{
//    int i;
//    double sumlog = 0.0, sum = 0.0;
//    double C[10] = { -6.089, -17.164, -34.054, -5.914, -24.721, -14.986, -24.100, -10.708, -26.662, -22.179 };

//    /* objective function */
//    for (i = 0; i < 10; i++)
//        sumlog += x[i];
//    for (i = 0; i < 10; i++)
//        sum += x[i] * (C[i] + log (x[i] / sumlog));
//    f[0] = sum;
//}

///// Implementation of the constraint function.
//void cec2006_g14::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
//{
//    /* constraints h=0 */
//    c[0] = x[0] + 2.0 * x[1] + 2.0 * x[2] + x[5] + x[9] - 2.0;
//    c[1] = x[3] + 2.0 * x[4] + x[5] + x[6] - 1.0;
//    c[2] = x[2] + x[6] + x[7] + 2.0 * x[8] + x[9] - 1.0;
//}

//std::string cec2006_g14::get_name() const
//{
//    return "CEC2006 - g14";
//}

//// -------------------------------------------

///// Default constructor
//cec2006_g15::cec2006_g15():base(3,0,1,2,0)
//{
//}

///// Clone method.
//base_ptr cec2006_g15::clone() const
//{
//    return base_ptr(new cec2006_g15(*this));
//}

///// Implementation of the objective function.
//void cec2006_g15::objfun_impl(fitness_vector &f, const decision_vector &x) const
//{
//    /* objective function */
//    f[0] = 1000.0 - pow(x[0], 2.0) - 2.0 * x[1] * x[1] - x[2] * x[2] - x[0] * x[1] - x[0] * x[2];
//}

///// Implementation of the constraint function.
//void cec2006_g15::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
//{
//    /* constraints h=0 */
//    h[0] = pow(x[0], 2.0) + pow(x[1], 2.0) + pow(x[2], 2.0) - 25.0;
//    h[1] = 8.0 * x[0] + 14.0 * x[1] + 7.0 * x[2] - 56.0;
//}

//std::string cec2006_g15::get_name() const
//{
//    return "CEC2006 - g15";
//}

//// -------------------------------------------

///// Default constructor
//cec2006_g16::cec2006_g16():base(5,0,1,38,38)
//{
//}

///// Clone method.
//base_ptr cec2006_g16::clone() const
//{
//    return base_ptr(new cec2006_g16(*this));
//}

///// Implementation of the objective function.
//void cec2006_g16::objfun_impl(fitness_vector &f, const decision_vector &x) const
//{
//    double x1, x2, x3, x4, x5;
//    double C[17], Y[17];

//    x1 = x[0];
//    x2 = x[1];
//    x3 = x[2];
//    x4 = x[3];
//    x5 = x[4];

//    Y[0] = x2 + x3 + 41.6;
//    C[0] = 0.024 * x4 - 4.62;
//    Y[1] = (12.5 / C[0]) + 12.0;
//    C[1] = 0.0003535 * pow (x1, 2.0) + 0.5311 * x1 + 0.08705 * Y[1] * x1;
//    C[2] = 0.052 * x1 + 78.0 + 0.002377 * Y[1] * x1;
//    Y[2] = C[1] / C[2];
//    Y[3] = 19.0 * Y[2];
//    C[3] = 0.04782 * (x1 - Y[2]) + ((0.1956 * pow (x1 - Y[2], 2.0)) / x2) + 0.6376 * Y[3] + 1.594 * Y[2];
//    C[4] = 100 * x2;
//    C[5] = x1 - Y[2] - Y[3];
//    C[6] = 0.950 - (C[3] / C[4]);
//    Y[4] = C[5] * C[6];
//    Y[5] = x1 - Y[4] - Y[3] - Y[2];
//    C[7] = (Y[4] + Y[3]) * 0.995;
//    Y[6] = C[7] / Y[0];
//    Y[7] = C[7] / 3798.0;
//    C[8] = Y[6] - (0.0663 * Y[6] / Y[7]) - 0.3153;
//    Y[8] = (96.82 / C[8]) + 0.321 * Y[0];
//    Y[9] = 1.29 * Y[4] + 1.258 * Y[3] + 2.29 * Y[2] + 1.71 * Y[5];
//    Y[10] = 1.71 * x1 - 0.452 * Y[3] + 0.580 * Y[2];
//    C[9] = 12.3 / 752.3;
//    C[10] = 1.75 * Y[1] * 0.995 * x1;
//    C[11] = 0.995 * Y[9] + 1998.0;
//    Y[11] = C[9] * x1 + (C[10] / C[11]);
//    Y[12] = C[11] - 1.75 * Y[1];
//    Y[13] = 3623.0 + 64.4 * x2 + 58.4 * x3 + (146312.0 / (Y[8] + x5));
//    C[12] = 0.995 * Y[9] + 60.8 * x2 + 48 * x4 - 0.1121 * Y[13] - 5095.0;
//    Y[14] = Y[12] / C[12];
//    Y[15] = 148000.0 - 331000.0 * Y[14] + 40.0 * Y[12] - 61.0 * Y[14] * Y[12];
//    C[13] = 2324 * Y[9] - 28740000 * Y[1];
//    Y[16] = 14130000 - 1328.0 * Y[9] - 531.0 * Y[10] + (C[13] / C[11]);
//    C[14] = (Y[12] / Y[14]) - (Y[12] / 0.52);
//    C[15] = 1.104 - 0.72 * Y[14];
//    C[16] = Y[8] + x5;

//    /* objective function */
//    f[0] = 0.0000005843 * Y[16] - 0.000117 * Y[13] - 0.1365 - 0.00002358 * Y[12] - 0.000001502 * Y[15] - 0.0321 * Y[11] - 0.004324 * Y[4] - 0.0001 * (C[14] / C[15]) - 37.48 * (Y[1] / C[11]);
//    f[0] = -f[0]; /* Max-->Min, Modified by Jane,Nov 22 2005 */
//}

///// Implementation of the constraint function.
//void cec2006_g16::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
//{
//    double Y[17];

//    Y[0] = x2 + x3 + 41.6;
//    C[0] = 0.024 * x4 - 4.62;
//    Y[1] = (12.5 / C[0]) + 12.0;
//    C[1] = 0.0003535 * pow (x1, 2.0) + 0.5311 * x1 + 0.08705 * Y[1] * x1;
//    C[2] = 0.052 * x1 + 78.0 + 0.002377 * Y[1] * x1;
//    Y[2] = C[1] / C[2];
//    Y[3] = 19.0 * Y[2];
//    C[3] = 0.04782 * (x1 - Y[2]) + ((0.1956 * pow (x1 - Y[2], 2.0)) / x2) + 0.6376 * Y[3] + 1.594 * Y[2];
//    C[4] = 100 * x2;
//    C[5] = x1 - Y[2] - Y[3];
//    C[6] = 0.950 - (C[3] / C[4]);
//    Y[4] = C[5] * C[6];
//    Y[5] = x1 - Y[4] - Y[3] - Y[2];
//    C[7] = (Y[4] + Y[3]) * 0.995;
//    Y[6] = C[7] / Y[0];
//    Y[7] = C[7] / 3798.0;
//    C[8] = Y[6] - (0.0663 * Y[6] / Y[7]) - 0.3153;
//    Y[8] = (96.82 / C[8]) + 0.321 * Y[0];
//    Y[9] = 1.29 * Y[4] + 1.258 * Y[3] + 2.29 * Y[2] + 1.71 * Y[5];
//    Y[10] = 1.71 * x1 - 0.452 * Y[3] + 0.580 * Y[2];
//    C[9] = 12.3 / 752.3;
//    C[10] = 1.75 * Y[1] * 0.995 * x1;
//    C[11] = 0.995 * Y[9] + 1998.0;
//    Y[11] = C[9] * x1 + (C[10] / C[11]);
//    Y[12] = C[11] - 1.75 * Y[1];
//    Y[13] = 3623.0 + 64.4 * x2 + 58.4 * x3 + (146312.0 / (Y[8] + x5));
//    C[12] = 0.995 * Y[9] + 60.8 * x2 + 48 * x4 - 0.1121 * Y[13] - 5095.0;
//    Y[14] = Y[12] / C[12];
//    Y[15] = 148000.0 - 331000.0 * Y[14] + 40.0 * Y[12] - 61.0 * Y[14] * Y[12];
//    C[13] = 2324 * Y[9] - 28740000 * Y[1];
//    Y[16] = 14130000 - 1328.0 * Y[9] - 531.0 * Y[10] + (C[13] / C[11]);
//    C[14] = (Y[12] / Y[14]) - (Y[12] / 0.52);
//    C[15] = 1.104 - 0.72 * Y[14];
//    C[16] = Y[8] + x5;

//    /* constraints g(x) <= 0 */
//    g[0] = -Y[3] + (0.28 / 0.72) * Y[4];
//    g[1] = -1.5 * x2 + x3;
//    g[2] = -21.0 + 3496.0 * (Y[1] / C[11]);
//    g[3] = -(62212.0 / C[16]) + 110.6 + Y[0];
//    g[4] = 213.1 - Y[0];
//    g[5] = Y[0] - 405.23;
//    g[6] = 17.505 - Y[1];
//    g[7] = Y[1] - 1053.6667;
//    g[8] = 11.275 - Y[2];
//    g[9] = Y[2] - 35.03;
//    g[10] = 214.228 - Y[3];
//    g[11] = Y[3] - 665.585;
//    g[12] = 7.458 - Y[4];
//    g[13] = Y[4] - 584.463;
//    g[14] = 0.961 - Y[5];
//    g[15] = Y[5] - 265.916;
//    g[16] = 1.612 - Y[6];
//    g[17] = Y[6] - 7.046;
//    g[18] = 0.146 - Y[7];
//    g[19] = Y[7] - 0.222;
//    g[20] = 107.99 - Y[8];
//    g[21] = Y[8] - 273.366;
//    g[22] = 922.693 - Y[9];
//    g[23] = Y[9] - 1286.105;
//    g[24] = 926.832 - Y[10];
//    g[25] = Y[10] - 1444.046;
//    g[26] = 18.766 - Y[11];
//    g[27] = Y[11] - 537.141;
//    g[28] = 1072.163 - Y[12];
//    g[29] = Y[12] - 3247.039;
//    g[30] = 8961.448 - Y[13];
//    g[31] = Y[13] - 26844.086;
//    g[32] = 0.063 - Y[14];
//    g[33] = Y[14] - 0.386;
//    g[34] = 71084.33 - Y[15];
//    g[35] = Y[15] - 140000.0;
//    g[36] = 2802713.0 - Y[16];
//    g[37] = Y[16] - 12146108.0;
//}

//std::string cec2006_g16::get_name() const
//{
//    return "CEC2006 - g16";
//}

// -------------------------------------------

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::cec2006);
