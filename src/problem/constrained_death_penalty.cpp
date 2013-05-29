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
#include <algorithm>

#include "../exceptions.h"
#include "../types.h"
#include "../population.h"
#include "base.h"
#include "constrained_death_penalty.h"

namespace pagmo { namespace problem {

/**
 * Constructor using initial constrained problem
 *
 * @param[in] problem base::problem to be modified to use a death penalty
 * as constraints handling technique.
 * @param[in] death_penalty_methos int to be modified to use a simple death penalty
 * if defined with 0 and a Kuri death penalty with a 1.
 *
 * @see Coello Coello, C. A. (2002). Theoretical and numerical constraint-handling techniques used with evolutionary algorithms: a survey of the state of the art. Computer methods in applied mechanics and engineering, 191(11), 1245-1287.
 */
constrained_death_penalty::constrained_death_penalty(const base &problem, const int death_penalty_method):
    base((int)problem.get_dimension(),
         problem.get_i_dimension(),
         problem.get_f_dimension(),
         0,
         0,
         0.),
    m_original_problem(problem.clone())
{
    if(m_original_problem->get_c_dimension() <= 0){
        pagmo_throw(value_error,"The original problem has no constraints.");
    }

    this->set_bounds(m_original_problem->get_lb(),m_original_problem->get_ub());

    // sets death penalty method
    this->set_death_penalty_method(death_penalty_method);
}

/// Copy Constructor. Performs a deep copy
constrained_death_penalty::constrained_death_penalty(const constrained_death_penalty &prob):
    base((int)prob.get_dimension(),
         prob.get_i_dimension(),
         prob.get_f_dimension(),
         prob.get_c_dimension(),
         prob.get_ic_dimension(),
         prob.get_c_tol()),
    m_original_problem(prob.m_original_problem->clone()),
    m_death_penalty_method(prob.m_death_penalty_method)
{
    this->set_bounds(m_original_problem->get_lb(),m_original_problem->get_ub());
}

/// Clone method.
base_ptr constrained_death_penalty::clone() const
{
    return base_ptr(new constrained_death_penalty(*this));
}

/// Sets the death penalty method
void constrained_death_penalty::set_death_penalty_method(const int death_penalty_method)
{
    if( death_penalty_method > 1 || death_penalty_method < 0) {
        pagmo_throw(value_error, "the death penalty method must be set to 0 for simple death or to 1 for Kuri death.");
    } else {
        this->m_death_penalty_method = death_penalty_method;
    }
}

/// Implementation of the objective function.
/// (Wraps over the original implementation)
void constrained_death_penalty::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
    constraint_vector c(m_original_problem->get_c_dimension(),0);
    m_original_problem->compute_constraints(c,x);

    if(m_original_problem->feasibility_c(c)) {
        m_original_problem->objfun(f, x);
    } else {
        double high_value = boost::numeric::bounds<double>::highest();

        switch(m_death_penalty_method)
        {
        case 0:
        {
            std::fill(f.begin(),f.end(),high_value);
            break;
        }
        case 1:
        {
            int number_of_constraints = c.size();
            int number_of_satisfied_constraints = 0;

            // computes the number of satisfied constraints
            for(c_size_type i=0; i<number_of_constraints; i++){
                if(m_original_problem->test_constraint(c,i))
                    number_of_satisfied_constraints += 1;
            }

            // sets the Kuri penalization
            double penalization = high_value * (1. - (double)number_of_satisfied_constraints / (double)number_of_constraints);

            std::fill(f.begin(),f.end(),penalization);
            break;
        }
        default:
            pagmo_throw(value_error, "Error: There are only 2 methods for the death penalty!");
            break;
        }
    }
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the type of constraint handling
 */
std::string constrained_death_penalty::human_readable_extra() const
{
    std::ostringstream oss;
    oss << m_original_problem->human_readable_extra() << std::endl;
    oss << "\n\tConstraints handled with death penalty." << std::endl;

    return oss.str();
}

std::string constrained_death_penalty::get_name() const
{
    return m_original_problem->get_name() + " [constrained_death_penalty, method_" +
            boost::lexical_cast<std::string>(m_death_penalty_method) + "]";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::constrained_death_penalty);

