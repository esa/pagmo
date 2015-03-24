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
 *   the Free Software Foundation; either version 3 of the License, or       *
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

#ifndef PAGMO_PROBLEM_PYTHON_BASE_H
#define PAGMO_PROBLEM_PYTHON_BASE_H

#include <boost/numeric/conversion/cast.hpp>
#include <boost/python/class.hpp>
#include <string>

#include "../../src/config.h"
#include "../../src/exceptions.h"
#include "../../src/problem/base.h"
#include "../../src/serialization.h"
#include "../../src/types.h"
#include "../utils.h"

// Forward declarations.
namespace pagmo { namespace problem {

class python_base;

}}

namespace boost { namespace serialization {

template <class Archive>
void save_construct_data(Archive &, const pagmo::problem::python_base *, const unsigned int);

template <class Archive>
inline void load_construct_data(Archive &, pagmo::problem::python_base *, const unsigned int);

}}

namespace pagmo { namespace problem {

class __PAGMO_VISIBLE python_base: public base, public boost::python::wrapper<base>
{
	public:
		explicit python_base(int n, int ni = 0, int nf = 1, int nc = 0, int nic = 0, const double &c_tol = 0):
			base(n,ni,nf,nc,nic,c_tol), boost::python::wrapper<base>() {}
		explicit python_base(int n, int ni, int nf, int nc, int nic, const std::vector<double> &c_tol):
			base(n,ni,nf,nc,nic,c_tol), boost::python::wrapper<base>() {}
		explicit python_base(const decision_vector &lb, const decision_vector &ub, int ni = 0, int nf = 1, int nc = 0, int nic = 0, const double &c_tol = 0):
			base(lb,ub,ni,nf,nc,nic,c_tol), boost::python::wrapper<base>() {}
		base_ptr clone() const
		{
			base_ptr retval = this->get_override("__get_deepcopy__")();
			if (!retval) {
				pagmo_throw(std::runtime_error,"problem's __get_deepcopy__() method returns a NULL pointer, please check the implementation");
			}
			return retval;
		}
		std::string get_name() const
		{
			if (boost::python::override f = this->get_override("get_name")) {
			#if BOOST_WORKAROUND(BOOST_MSVC, <= 1700)
				return boost::python::call<std::string>(this->get_override("get_name").ptr());
			#else
				return f();				
			#endif
			}
			return base::get_name();
		}
		std::string default_get_name() const
		{
			return this->base::get_name();
		}
		std::string human_readable_extra() const
		{
			if (boost::python::override f = this->get_override("human_readable_extra")) {
			#if BOOST_WORKAROUND(BOOST_MSVC, <= 1700)
				return boost::python::call<std::string>(this->get_override("human_readable_extra").ptr());
			#else
				return f();				
			#endif

			}
			return base::human_readable_extra();
		}
		std::string default_human_readable_extra() const
		{
			return this->base::human_readable_extra();
		}
		fitness_vector py_objfun(const decision_vector &x) const
		{
			if (boost::python::override f = this->get_override("_objfun_impl")) {
				return f(x);
			}
			pagmo_throw(not_implemented_error,"objective function has not been implemented");
		}
		std::string get_typename() const
		{
			if (boost::python::override f = this->get_override("_get_typename")) {
			#if BOOST_WORKAROUND(BOOST_MSVC, <= 1700)
				return boost::python::call<std::string>(this->get_override("_get_typename").ptr());
			#else
				return f();				
			#endif
			}
			pagmo_throw(not_implemented_error,"the '_get_typename' method has not been implemented");
		}
		bool py_equality_operator_extra(const base &p) const
		{
			if (boost::python::override f = this->get_override("_equality_operator_extra")) {
				return f(p);
			}
			return base::equality_operator_extra(p);
		}
		constraint_vector py_compute_constraints_impl(const decision_vector &x) const
		{
			boost::python::override f = this->get_override("_compute_constraints_impl");
			pagmo_assert(f);
			return f(x);
        }
        bool py_compare_fitness_impl(const fitness_vector &f0, const fitness_vector &f1) const
        {
            boost::python::override f = this->get_override("_compare_fitness_impl");
            pagmo_assert(f);
            return f(f0, f1);
        }
        bool py_compare_constraints_impl(const constraint_vector &c0, const constraint_vector &c1) const
        {
            boost::python::override f = this->get_override("_compare_constraints_impl");
            pagmo_assert(f);
            return f(c0, c1);
        }
        bool py_compare_fc_impl(const fitness_vector &f0, const constraint_vector &c0, const fitness_vector &f1, const constraint_vector &c1) const
        {
            boost::python::override f = this->get_override("_compare_fc_impl");
            pagmo_assert(f);
            return f(f0, c0, f1, c1);
        }

	protected:
		void objfun_impl(fitness_vector &f, const decision_vector &x) const
		{
			f = py_objfun(x);
		}
		bool equality_operator_extra(const base &p) const
		{
			// NOTE: here the dynamic cast is safe because in base equality we already checked the C++ type.
			if (get_typename() != dynamic_cast<const python_base &>(p).get_typename()) {
				return false;
			}
			return py_equality_operator_extra(p);
		}
		void compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
		{
			if (this->get_override("_compute_constraints_impl")) {
				// If the function is overridden, use it.
				c = py_compute_constraints_impl(x);
			} else {
				// Otherwise, avoid memory allocation by calling the base function directly.
				base::compute_constraints_impl(c,x);
			}
		}
        bool compare_fitness_impl(const fitness_vector &f0, const fitness_vector &f1) const
        {
            if(this->get_override("_compare_fitness_impl")) {
                // if the function is overidden, use it
                return py_compare_fitness_impl(f0, f1);
            } else {
                // else, the base function is called directly
                return problem::base::compare_fitness_impl(f0, f1);
            }
        }
        bool compare_constraints_impl(const constraint_vector &c0, const constraint_vector &c1) const
        {
            if(this->get_override("_compare_constraints_impl")) {
                // if the function is overidden, use it
                return py_compare_constraints_impl(c0, c1);
            } else {
                // else, the base function is called directly
                return problem::base::compare_constraints_impl(c0, c1);
            }
        }
        bool compare_fc_impl(const fitness_vector &f0, const constraint_vector &c0, const fitness_vector &f1, const constraint_vector &c1) const
        {
            if(this->get_override("_compare_fc_impl")) {
                // if the function is overidden, use it
                return py_compare_fc_impl(f0, c0, f1, c1);
            } else {
                // else, the base function is called directly
                return problem::base::compare_fc_impl(f0, c0, f1, c1);
            }
        }

	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & boost::serialization::base_object<boost::python::wrapper<base> >(*this);
		}
};

} }

namespace boost { namespace serialization {

template <class Archive>
inline void save_construct_data(Archive &ar, const pagmo::problem::python_base *prob, const unsigned int)
{
	// Save data required to construct instance.
	int dimension = boost::numeric_cast<int>(prob->get_dimension());
	ar << dimension;
}

template <class Archive>
inline void load_construct_data(Archive &ar, pagmo::problem::python_base *prob, const unsigned int)
{
	// Retrieve data from archive required to construct new instance.
	int dimension;
	ar >> dimension;
	// Invoke inplace constructor to initialize instance of the island.
	::new(prob)pagmo::problem::python_base(dimension);
}

}} //namespaces

BOOST_CLASS_EXPORT(pagmo::problem::python_base)

#endif
