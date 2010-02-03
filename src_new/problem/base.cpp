/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
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

// 30/01/10 Created by Francesco Biscani.

#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/ref.hpp>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

#include "../atomic_counters/atomic_counters.h"
#include "../exceptions.h"
#include "../types.h"
#include "base.h"

namespace pagmo
{
namespace problem {

// Initialisation of static objective function calls counter.
atomic_counter_size_t base::m_objfun_counter(0);

/// Constructor from global, integer and fitness dimensions.
/**
 * Initialises global dimension to n, integer dimension to ni and fitness dimension to nf. n and nf must be positive and ni must be in the [0,n] range.
 * Lower and upper bounds are set to 0 and 1 respectively.
 */
base::base(int n, int ni, int nf):m_decision_vector_cache(cache_capacity),m_fitness_vector_cache(cache_capacity)
{
	if (n <= 0 || ni < 0 || nf <= 0 || ni > n) {
		pagmo_throw(value_error,"invalid dimension(s)");
	}
	m_i_dimension = boost::numeric_cast<size_type>(ni);
	m_f_dimension = boost::numeric_cast<f_size_type>(nf);
	const size_type size = boost::numeric_cast<size_type>(n);
	m_lb.resize(size);
	m_ub.resize(size);
	std::fill(m_lb.begin(),m_lb.end(),0);
	std::fill(m_ub.begin(),m_ub.end(),1);
	// Resize properly temporary fitness storage.
	m_tmp_f1.resize(m_f_dimension);
	m_tmp_f2.resize(m_f_dimension);
}

/// Constructor from upper/lower bounds and integer and fitness dimensions.
/**
 * Will fail if ni is negative or greater than lb.size(), if nf is not positive, if the sizes of the lower/upper bounds are zero or not identical, or
 * any lower bound is greater than the corresponding upper bound.
 */
base::base(const decision_vector &lb, const decision_vector &ub, int ni, int nf):m_decision_vector_cache(cache_capacity),m_fitness_vector_cache(cache_capacity)
{
	if (ni < 0 || nf <= 0 || boost::numeric_cast<size_type>(ni) > lb.size()) {
		pagmo_throw(value_error,"invalid dimension(s)");
	}
	m_i_dimension = boost::numeric_cast<size_type>(ni);
	m_f_dimension = boost::numeric_cast<f_size_type>(nf);
	construct_from_iterators(lb.begin(),lb.end(),ub.begin(),ub.end());
	// Resize properly temporary fitness storage.
	m_tmp_f1.resize(m_f_dimension);
	m_tmp_f2.resize(m_f_dimension);
}

/// Trivial destructor.
base::~base() {std::cout << m_decision_vector_cache.size() << '\n'; std::cout << m_fitness_vector_cache.size() << '\n';}

/// Lower bounds getter.
const decision_vector &base::get_lb() const
{
	return m_lb;
}

/// Upper bounds getter.
const decision_vector &base::get_ub() const
{
	return m_ub;
}

/// Bounds setter from pagmo::decision_vector.
/**
 * Set lower/upper bounds to lb/ub. Will fail if lb and ub sizes do not match, if their sizes are different
 * from the global size of the problem or if at least one lower bound is greater than the corresponding upper bound.
 */
void base::set_bounds(const decision_vector &lb, const decision_vector &ub)
{
	if (lb.size() != ub.size() || lb.size() != m_lb.size()) {
		pagmo_throw(value_error,"invalid or inconsistent bounds dimensions in set_bounds()");
	}
	verify_bounds(lb.begin(),lb.end(),ub.begin(),ub.end());
	m_lb = lb;
	m_ub = ub;
}

/// Set lower bounds from pagmo::decision_vector.
/**
 * Will fail if lb's size is different from the global size or if at least one lower bound is greater than the corresponding upper bound.
 */
void base::set_lb(const decision_vector &lb)
{
	if (lb.size() != m_lb.size()) {
		pagmo_throw(value_error,"invalid bounds dimension in set_lb()");
	}
	verify_bounds(lb.begin(),lb.end(),m_ub.begin(),m_ub.end());
	m_lb = lb;
}

/// Set specific lower bound to value.
/**
 * Will fail if n overflows global dimension or if value is greater than the corresponding upper bound.
 */
void base::set_lb(int n, const double &value)
{
	const size_type i = boost::numeric_cast<size_type>(n);
	if (i >= m_lb.size() || m_ub[i] < value) {
		pagmo_throw(value_error,"invalid index and/or value for lower bound");
	}
	m_lb[i] = value;
}

/// Set all lower bounds to value.
/**
 * Will fail if value is greater than at least one upper bound.
 */
void base::set_lb(const double &value)
{
	for (size_type i = 0; i < m_lb.size(); ++i) {
		if (m_ub[i] < value) {
			pagmo_throw(value_error,"invalid value for lower bound");
		}
	}
	std::fill(m_lb.begin(),m_lb.end(),value);
}

/// Set upper bounds from pagmo::decision_vector.
/**
 * Will fail if ub's size is different from the global size or if at least one upper bound is less than the corresponding lower bound.
 */
void base::set_ub(const decision_vector &ub)
{
	if (ub.size() != m_lb.size()) {
		pagmo_throw(value_error,"invalid bounds dimension in set_ub()");
	}
	verify_bounds(m_lb.begin(),m_lb.end(),ub.begin(),ub.end());
	m_ub = ub;
}

/// Set specific upper bound to value.
/**
 * Will fail if n overflows global dimension or if value is less than the corresponding lower bound.
 */
void base::set_ub(int n, const double &value)
{
	const size_type i = boost::numeric_cast<size_type>(n);
	if (i >= m_lb.size() || m_lb[i] > value) {
		pagmo_throw(value_error,"invalid index and/or value for upper bound");
	}
	m_ub[i] = value;
}

/// Set all upper bounds to value.
/**
 * Will fail if value is less than at least one lower bound.
 */
void base::set_ub(const double &value)
{
	for (size_type i = 0; i < m_lb.size(); ++i) {
		if (m_lb[i] > value) {
			pagmo_throw(value_error,"invalid value for upper bound");
		}
	}
	std::fill(m_ub.begin(),m_ub.end(),value);
}

/// Return global dimension.
base::size_type base::get_dimension() const
{
	return m_lb.size();
}

/// Return integer dimension.
base::size_type base::get_i_dimension() const
{
	return m_i_dimension;
}

/// Return fitness dimension.
base::f_size_type base::get_f_dimension() const
{
	return m_f_dimension;
}

/// Return fitness of pagmo::decision_vector.
/**
 * Equivalent to:
 @verbatim
 fitness_vector f(get_f_dimension());
 objfun(f,x);
 return f;
 @endverbatim
 */
fitness_vector base::objfun(const decision_vector &x) const
{
	fitness_vector f(m_f_dimension);
	objfun(f,x);
	return f;
}

/// Write fitness of pagmo::decision_vector into pagmo::fitness_vector.
/**
 * Will call objfun_impl() internally. Will fail if f's and x's size are different from the fitness and global dimension respectively
 * or if the decision vector is outside the bounds of the problem.
 *
 * The implementation internally uses a caching mechanism, so that recently-computed quantities are remembered and re-used when appropriate.
 */
void base::objfun(fitness_vector &f, const decision_vector &x) const
{
	// Some checks on the input values.
	if (f.size() != m_f_dimension || x.size() != m_lb.size()) {
		pagmo_throw(value_error,"vector size(s) mismatch(es) when calling objective function");
	}
	for (size_type i = 0; i < m_lb.size(); ++i) {
		if (x[i] < m_lb[i] || x[i] > m_ub[i]) {
			pagmo_throw(value_error,"decision vector is not within problem's bounds");
		}
	}
	// Look into the cache.
	typedef decision_vector_cache_type::iterator x_iterator;
	typedef fitness_vector_cache_type::iterator f_iterator;
	const x_iterator x_it = std::find(m_decision_vector_cache.begin(),m_decision_vector_cache.end(),x);
	if (x_it == m_decision_vector_cache.end()) {
std::cout << "cache miss\n";
		// Fitness is not into memory. Calculate it.
		objfun_impl(f,x);
		// Store the decision vector and the newly-calculated fitness in the front of the buffers.
		m_decision_vector_cache.push_front(x);
		m_fitness_vector_cache.push_front(f);
	} else {
std::cout << "cache hit\n";
		// Compute the corresponding iterator in the fitness vector cache.
		f_iterator f_it = m_fitness_vector_cache.begin();
		std::advance(f_it,std::distance(m_decision_vector_cache.begin(),x_it));
		// Assign to the fitness vector the value in the cache.
		f = *f_it;
		// Move the content of the current positions to the front of the buffers, and shift everything else down
		// by one position.
		x_iterator tmp_x_it = m_decision_vector_cache.begin();
		f_iterator tmp_f_it = m_fitness_vector_cache.begin();
		while (x_it != tmp_x_it) {
			x_it->swap(*tmp_x_it);
			f_it->swap(*tmp_f_it);
			++tmp_x_it;
			++tmp_f_it;
		}
		pagmo_assert(tmp_f_it == f_it);
	}
	// Actually do the increment only if we have fast incrementing capabilities in m_objfun_counter.
	if (m_objfun_counter.is_increment_fast) {
		++m_objfun_counter;
	}
}

/// Return human readable representation of the problem.
/**
 * Will return a formatted string containing:
 * - problem type (in mangled C++ form),
 * - dimensions,
 * - lower and upper bounds.
 *
 * The output of human_readable_extra() will be appended at the end of the string.
 */
std::string base::human_readable() const
{
	std::ostringstream s;
	s << "Problem type: " << typeid(*this).name() << '\n';
	const size_type size = get_dimension();
	s << "Global dimension:\t" << size << '\n';
	s << "Integer dimension:\t" << m_i_dimension << '\n';
	s << "Fitness dimension:\t" << m_f_dimension << '\n';
	s << "Lower bounds: ";
	s << m_lb << '\n';
	s << "Upper bounds: ";
	s << m_ub << '\n';
	s << human_readable_extra();
	return s.str();
}

/// Extra information in human readable format.
/**
 * Default implementation returns an empty string.
 */
std::string base::human_readable_extra() const
{
	return std::string();
}

/// Equality operator.
/**
 * The following conditions will be tested, in order:
 * - problems are of the same type,
 * - problems have the same dimension,
 * - lower and upper bounds are equal,
 * - return value of equality_operator_extra().
 *
 * If any of the conditions above is false, then the return value will also be false. Otherwise return value will be true.
 */
bool base::operator==(const base &p) const
{
	const size_type size = get_dimension();
	if (typeid(*this) != typeid(p) || size != p.get_dimension() || m_i_dimension != p.m_i_dimension || m_f_dimension != p.m_f_dimension) {
		return false;
	}
	for (size_t i = 0; i < size; ++i) {
		if (m_lb[i] != p.m_lb[i] || m_ub[i] != p.m_ub[i]) {
			return false;
		}
	}
	if (!equality_operator_extra(p)) {
		return false;
	}
	return true;
}

/// Compare decision vectors.
/**
 * This functions returns true if x1 is a better decision_vector than x2, false otherwise. This function will compute the
 * fitness vectors associated to x1 and x2 and will feed them to compare_impl(), whose result will be returned.
 */
bool base::compare(const decision_vector &x1, const decision_vector &x2) const
{
	// Make sure the size of the tmp fitness vectors are suitable.
	pagmo_assert(m_tmp_f1.size() == m_f_dimension && m_tmp_f2.size() == m_f_dimension);
	// Store fitnesses into temporary space.
	objfun(m_tmp_f1,x1);
	objfun(m_tmp_f2,x2);
	// Call the comparison implementation.
	return compare_impl(m_tmp_f1,m_tmp_f2);
}

/// Compare fitness vectors.
/**
 * Default implementation will compute the summations f1 and f2 of all elements of the input fitness vectors and will return f1 < f2.
 */
bool base::compare_impl(const fitness_vector &v_f1, const fitness_vector &v_f2) const
{
	typedef fitness_vector::value_type fitness_type;
	const fitness_type init = 0;
	const fitness_type f1 = std::accumulate(v_f1.begin(),v_f1.end(),init);
	const fitness_type f2 = std::accumulate(v_f2.begin(),v_f2.end(),init);
	return (f1 < f2);
}

/// Extra requirements for equality.
/**
 * Additional problem-specific equality testing. Default implementation returns true.
 */
bool base::equality_operator_extra(const base &) const
{
	return true;
}

/// Inequality operator.
/**
 * Equivalent to the negation of equality operator.
 */
bool base::operator!=(const base &p) const
{
	return !(*this == p);
}

/// Overload stream operator for problem::base.
/**
 * Equivalent to printing base::human_readable() to stream.
 */
std::ostream &operator<<(std::ostream &s, const base &p)
{
	s << p.human_readable();
	return s;
}

/// Return the total number of calls to the objective function.
/**
 * The number is a static global variable that gets incremented each time base::objfun() is called.
 */
std::size_t objfun_calls()
{
	if (!base::m_objfun_counter.is_increment_fast) {
		pagmo_throw(not_implemented_error,"fast atomic counters are not available in this version of PaGMO");
	}
	return (base::m_objfun_counter).get_value();
}

/// Reset to zero the total number of calls to the objective function.
void reset_objfun_calls()
{
	base::m_objfun_counter = atomic_counter_size_t();
}

}
}
