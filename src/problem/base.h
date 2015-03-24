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

// 30/01/10 Created by Francesco Biscani.

#ifndef PAGMO_PROBLEM_BASE_H
#define PAGMO_PROBLEM_BASE_H

// This define was added for a MSVC compilation fix and later removed as it creates a potential
// problem when testing pagmo in debug
// #define BOOST_CB_DISABLE_DEBUG 

#include <algorithm>
#include <boost/circular_buffer.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/shared_ptr.hpp>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <string>

#include "../config.h"
#include "../exceptions.h"
#include "../serialization.h"
#include "../types.h"
//#include "base_meta.h"

namespace pagmo
{
// Fwd declarations.
class population;
class island;

/// Problem namespace.
/**
 * This namespace contains all the problems implemented in PaGMO.
 */
namespace problem {

/// Base problem class.
class base;

/// Alias for shared pointer to base problem.
typedef boost::shared_ptr<base> base_ptr;

/// Base problem class.
/**
 * \section Introduction
 * This class represents a box-bounded, multiobjective, mixed-integer, constrained optimisation problem defined by:
 * - a global dimension, i.e., the number of dimensions of the global search space,
 * - the dimension of the integral (or combinatorial) part of the problem,
 * - lower and upper bounds of the global search space,
 * - the total number of constraints,
 * - the number of inequality constraints (never exceeding the total number of constraints),
 * - a constraint computation function,
 * - an objective function that takes as input a decision vector and returns a vector of fitnesses,
 * - a fitness dimension, i.e., the length of the fitness vector returned by the objective function,
 * - a constraints tolerance vector.
 *
 * All dimensions are invariant in the life cycle of a problem object.
 *
 * The bounds of the problem are allowed to vary over the whole range of double-precision values for continuous optimisation,
 * while for combinatorial optimisation the bounds must be in the [-32767,32767] range (corresponding to the INT_MIN and INT_MAX
 * constants defined in the C++ standard "climits" header file). All bounds setting functions will make sure that the following conditions are
 * respected:
 * - lower bounds are not greater than upper bounds,
 * - the bounds of the integer part of the problem are integer and they are within the allowed range.
 *
 * If the first condition is not met, an error will be raised. If the second condition is not met, the bounds will be set to the extremes
 * of the allowed range and/or rounded to the nearest integer as necessary. After that, an error will be generated in order to alert the user.
 *
 * The constraint tolerance vector contains the tolerances for each constraints. This allows to implement constraints of type:
 * - equality constraints h_i(x) < tol_i
 * - inequality constraints g_k(x) < tol_k
 *
 * Note that if a single value is given to the problem constructor instead of a constraint tolerance vector, then a constraint tolerance vector of size the number of constraints filled with the given value is created.
 *
 * All problems implemented in PaGMO must derive from this base class and implement the following pure virtual methods:
 * - the clone() method, i.e., the polymorphic copy constructor,
 * - the objfun_impl() method, i.e., the implementation of the objective function, which computes the fitness vector associated to a decision vector.
 *
 * Additionally, the following virtual protected methods can be reimplemented in derived classes:
 * - get_name(), for specifying a string identifier for the problem type,
 * - human_readable_extra(), for providing extra output when printing the problem to stream,
 * - equality_operator_extra(), for providing additional criterions when testing for equality between two problems,
 * - compare_fitness_impl(), to reimplement the function that compares two fitness vectors (returning true if the first vector is strictly better
 *   than the second one, false otherwise),
 * - compute_constraints_impl(), to calculate the constraint vector associated to a decision vector,
 * - compare_constraints_impl(), to compare two constraint vectors,
 * - compare_fc_impl(), to perform a simultaneous fitness/constraint vector pairs comparison.
 *
 * Please note that while a problem is intended to provide methods for ranking decision and constraint vectors, such methods are not to be used
 * mandatorily by an algorithm: each algorithm can decide to use its own ranking schemes during an optimisation. The ranking methods provided
 * by the problem are always used instead during the migration of decision vectors from one island to the other.
 *
 * \section Caching
 * A caching mechanism is implemented to make sure the objective function is never evaluated twice on the very same chromosome
 *
 * \section Serialization
 * The problem classes are serialized for the purpose of transmitting their corresponding objects over a distributed environment, as being part of the population class.
 * Serializing a derived problem requires that the needed serialization libraries be declared in the header of the derived class.
 * Virtually all the derived problem classes need to have the following declared in their header files:
@verbatim
	friend class boost::serialization::access;
	template<class Archive>
@endverbatim
 * Each derived class must implement implicitly, in its header file, the serialize method, which must contain the pointer to the base class like:
@verbatim
	ar & boost::serialization::base_object<base>(*this);
@endverbatim
 * and the rest of the attributes simply as archive members:
@verbatim
	ar & attribute_name;
@endverbatim
 * In order to be able to identify the derived problem class when deserializing an object declared as a base_pointer, the derived problem class needs to be registered. In the case where your derived problem class has a default constructor, this is done by registering the class in the "pagmo/src/problems.h", in the REGISTER_PROBLEM_SERIALIZATIONS() routine, by adding:
@verbatim
ar.template register_type<problem::derived_problem>();
@endverbatim
 * In the case where the derived problem does not have a default constructor, the serialization functions "load_construct_data" and "save_construct_data" need to be overriden. This done by invoking the non-default constructor in-place (in the load_construct_data) to initilize the memory (Examples can be found in the implemented problems or in the Boost Serialization libraray unde "Non-default constructor" section).W
 * Notes:
 * - "const" attributes need to be cast as constants in the serialize method using const_cast
 * - attributes that that are not primitives, need be a serialized type as well
 * - pointers to primitives cannot be serialized (in this case one can split the serialize method into save/load methods and store the values, that the pointers refer to, into temporary variables which are serialized insted - see boost serialize documentation on the topic if needed)
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE base
{
		// Meta problems need to be able to access protected virtual functions
		friend class base_meta;
		// Underlying containers used for caching decision and fitness vectors.
		typedef boost::circular_buffer<decision_vector> decision_vector_cache_type;
		typedef boost::circular_buffer<fitness_vector> fitness_vector_cache_type;
		typedef boost::circular_buffer<constraint_vector> constraint_vector_cache_type;
	public:
		/// Capacity of the internal caches.
		static const std::size_t cache_capacity = 5;
		/// Problem's size type: the same as pagmo::decision_vector's size type.
		typedef decision_vector::size_type size_type;
		/// Fitness' size type: the same as pagmo::fitness_vector's size type.
		typedef fitness_vector::size_type f_size_type;
		/// Constraints' size type: the same as pagmo::constraint_vector's size type.
		typedef constraint_vector::size_type c_size_type;
		base(int, int = 0, int = 1, int = 0, int = 0, const double & = 0);
		base(int, int, int, int, int, const std::vector<double> &);
		base(const double &, const double &, int, int = 0, int = 1, int = 0, int = 0, const double & = 0);
		base(const decision_vector &, const decision_vector &, int = 0, int = 1, int = 0, int = 0, const double & = 0);
		/// Constructor from raw arrays, integer dimension, fitness dimension, global constraints dimension, inequality constraints dimension and constraints tolerance.
		/**
		 * Lower and upper bounds are initialised with the content of two arrays of size N.
		 * Construction will fail if at least one lower bound is greater than the corresponding upper bound, if N is zero,
		 * if integer dimension is either negative or greater than the global dimension, if fitness dimension is not positive,
		 * if constraints dimensions are negative, if inequality constraints dimension is greater than global constraints dimension, or
		 * if constraints tolerance is negative.
		 *
		 * @param[in] v1 lower bounds for the problem.
		 * @param[in] v2 upper bounds for the problem.
		 * @param[in] ni dimension of the combinatorial part of the problem.
		 * @param[in] nf dimension of the fitness vector of the problem.
		 * @param[in] nc global number of constraints.
		 * @param[in] nic number of inequality constraints.
		 * @param[in] c_tol constraints tolerance. Fills the tolerance vector of size nc with c_tol.
		 */
		template <std::size_t N>
		base(const double (&v1)[N], const double (&v2)[N], int ni = 0, int nf = 1, int nc = 0, int nic = 0, const double &c_tol = 0):
			m_i_dimension(boost::numeric_cast<size_type>(ni)),m_f_dimension(boost::numeric_cast<f_size_type>(nf)),
			m_c_dimension(boost::numeric_cast<c_size_type>(nc)),m_ic_dimension(boost::numeric_cast<c_size_type>(nic)),
			m_c_tol(nc,c_tol),
			m_decision_vector_cache_f(boost::numeric_cast<decision_vector_cache_type::size_type>(cache_capacity)),
			m_fitness_vector_cache(boost::numeric_cast<fitness_vector_cache_type::size_type>(cache_capacity)),
			m_decision_vector_cache_c(boost::numeric_cast<decision_vector_cache_type::size_type>(cache_capacity)),
			m_constraint_vector_cache(boost::numeric_cast<constraint_vector_cache_type::size_type>(cache_capacity))
		{
			if (c_tol < 0) {
				pagmo_throw(value_error,"constraints tolerance must be non-negative");
			}
			if (!m_f_dimension) {
				pagmo_throw(value_error,"fitness dimension must be strictly positive");
			}
			if (m_ic_dimension > m_c_dimension) {
				pagmo_throw(value_error,"inequality constraints dimension must not be greater than global constraints dimension");
			}
			construct_from_iterators(v1,v1 + N,v2,v2 + N);
			if (m_i_dimension > m_lb.size()) {
				pagmo_throw(value_error,"integer dimension must not be greater than global dimension");
			}
			// Resize properly temporary fitness and constraint storage.
			m_tmp_f1.resize(m_f_dimension);
			m_tmp_f2.resize(m_f_dimension);
			m_tmp_c1.resize(m_c_dimension);
			m_tmp_c2.resize(m_c_dimension);
			// Normalise bounds.
			normalise_bounds();
		}
		/// Constructor from iterators, integer dimension, fitness dimension, global constraints dimension, inequality constraints dimension and constraints tolerance.
		/**
		 * Lower bounds are initialised with the content in the range [start1,end1[, upper bounds with the content in the range [start2,end2[.
		 * Construction will fail if the ranges have different or null sizes, if at least one lower bound is greater than the corresponding upper bound,
		 * if integer dimension is either negative or greater than the global dimension, if fitness dimension is not positive,
		 * if constraints dimensions are negative, if inequality constraints dimension is greater than global constraints dimension or
		 * if constraints tolerance is negative.
		 *
		 * @param[in] start1 iterator to the beginning of the lower bounds sequence.
		 * @param[in] end1 iterator to the end of the lower bounds sequence.
		 * @param[in] start2 iterator to the beginning of the upper bounds sequence.
		 * @param[in] end2 iterator to the end of the upper bounds sequence.
		 * @param[in] ni dimension of the combinatorial part of the problem.
		 * @param[in] nf dimension of the fitness vector of the problem.
		 * @param[in] nc global number of constraints.
		 * @param[in] nic number of inequality constraints.
		 * @param[in] c_tol constraints tolerance. Fills the tolerance vector of size nc with c_tol..
		 */
		template <class Iterator1, class Iterator2>
		base(Iterator1 start1, Iterator1 end1, Iterator2 start2, Iterator2 end2, int ni = 0, int nf = 1, int nc = 0, int nic = 0, const double &c_tol = 0):
			m_i_dimension(boost::numeric_cast<size_type>(ni)),m_f_dimension(boost::numeric_cast<f_size_type>(nf)),
			m_c_dimension(boost::numeric_cast<c_size_type>(nc)),m_ic_dimension(boost::numeric_cast<c_size_type>(nic)),
			m_c_tol(nc,c_tol),
			m_decision_vector_cache_f(boost::numeric_cast<decision_vector_cache_type::size_type>(cache_capacity)),
			m_fitness_vector_cache(boost::numeric_cast<fitness_vector_cache_type::size_type>(cache_capacity)),
			m_decision_vector_cache_c(boost::numeric_cast<decision_vector_cache_type::size_type>(cache_capacity)),
			m_constraint_vector_cache(boost::numeric_cast<constraint_vector_cache_type::size_type>(cache_capacity))
		{
			if (c_tol < 0) {
				pagmo_throw(value_error,"constraints tolerance must be non-negative");
			}
			if (!m_f_dimension) {
				pagmo_throw(value_error,"fitness dimension must be strictly positive");
			}
			if (m_ic_dimension > m_c_dimension) {
				pagmo_throw(value_error,"inequality constraints dimension must not be greater than global constraints dimension");
			}
			construct_from_iterators(start1,end1,start2,end2);
			if (m_i_dimension > m_lb.size()) {
				pagmo_throw(value_error,"integer dimension must not be greater than global dimension");
			}
			// Properly resize temporary fitness and constraint storage.
			m_tmp_f1.resize(m_f_dimension);
			m_tmp_f2.resize(m_f_dimension);
			m_tmp_c1.resize(m_c_dimension);
			m_tmp_c2.resize(m_c_dimension);
			// Normalise bounds.
			normalise_bounds();
		}
		virtual ~base();
		/** @name Bounds setters/getters.
		* Methods used to manipulate problem bounds.
		*/
		//@{
		const decision_vector &get_lb() const;
		const decision_vector &get_ub() const;
		void set_bounds(const decision_vector &, const decision_vector &);
		/// Bounds setter from iterators.
		/**
		 * Set lower and upper bounds to the content of the ranges [start1,end1[ and [start2,end2[. Will fail if ranges sizes do not match, if ranges sizes are different
		 * from the global size of the problem or if at least one lower bound is greater than the corresponding upper bound.
		 *
		 * @param[in] start1 iterator to the beginning of the lower bounds sequence.
		 * @param[in] end1 iterator to the end of the lower bounds sequence.
		 * @param[in] start2 iterator to the beginning of the upper bounds sequence.
		 * @param[in] end2 iterator to the end of the upper bounds sequence.
		 */
		template <class Iterator1, class Iterator2>
		void set_bounds(Iterator1 start1, Iterator1 end1, Iterator2 start2, Iterator2 end2)
		{
			typedef typename std::iterator_traits<Iterator1>::difference_type d_type1;
			typedef typename std::iterator_traits<Iterator2>::difference_type d_type2;
			const d_type1 d1 = std::distance(start1,end1);
			const d_type2 d2 = std::distance(start2,end2);
			if (d1 != d2 || d1 != std::distance(m_lb.begin(),m_lb.end())) {
				pagmo_throw(value_error,"invalid or inconsistent bounds dimensions in set_bounds()");
			}
			verify_bounds(start1,end1,start2,end2);
			std::copy(start1,end1,m_lb.begin());
			std::copy(start2,end2,m_ub.begin());
			// Normalise bounds.
			normalise_bounds();
		}
		/// Bounds setter from raw arrays.
		/**
		 * Set lower and upper bounds to the content of the raw arrays v1 and v2. Will fail if N is different
		 * from the global size of the problem or if at least one lower bound is greater than the corresponding upper bound.
		 *
		 * @param[in] v1 lower bounds for the problem.
		 * @param[in] v2 upper bounds for the problem.
		 */
		template <std::size_t N>
		void set_bounds(const double (&v1)[N], const double (&v2)[N])
		{
			if (m_lb.size() != N) {
				pagmo_throw(value_error,"invalid bounds dimensions in set_bounds()");
			}
			verify_bounds(v1,v1 + N,v2,v2 + N);
			std::copy(v1,v1 + N,m_lb.begin());
			std::copy(v2,v2 + N,m_ub.begin());
			// Normalise bounds.
			normalise_bounds();
		}
		void set_bounds(const double &, const double &);
		void set_bounds(int, const double &, const double &);
		void set_lb(const decision_vector &);
		void set_lb(int, const double &);
		void set_lb(const double &);
		/// Lower bounds setter from iterators.
		/**
		 * Will fail if the iterator distance is different from global problem dimension or if at least one lower bound is greater than the
		 * corresponding upper bound.
		 *
		 * @param[in] start iterator to the beginning of the lower bounds sequence.
		 * @param[in] end iterator to the end of the lower bounds sequence.
		 */
		template <class Iterator>
		void set_lb(Iterator start, Iterator end)
		{
			if (std::distance(start,end) != std::distance(m_lb.begin(),m_lb.end())) {
				pagmo_throw(value_error,"invalid bounds dimension in set_lb()");
			}
			verify_bounds(start,end,m_ub.begin(),m_ub.end());
			std::copy(start,end,m_lb.begin());
			// Normalise bounds.
			normalise_bounds();
		}
		/// Lower bounds setter from raw array.
		/**
		 * Will fail if N is different from global problem dimension or if at least one lower bound is greater than the
		 * corresponding upper bound.
		 *
		 * @param[in] v lower bounds array.
		 */
		template <std::size_t N>
		void set_lb(const double (&v)[N])
		{
			if (N != m_lb.size()) {
				pagmo_throw(value_error,"invalid bounds dimension in set_lb()");
			}
			verify_bounds(v,v + N,m_ub.begin(),m_ub.end());
			std::copy(v,v + N,m_lb.begin());
			// Normalise bounds.
			normalise_bounds();
		}
		void set_ub(const decision_vector &);
		void set_ub(int, const double &);
		void set_ub(const double &);
		/// Upper bounds setter from iterators.
		/**
		 * Will fail if the iterator distance is different from global problem dimension or if at least one upper bound is less than the
		 * corresponding lower bound.
		 *
		 * @param[in] start iterator to the beginning of the upper bounds sequence.
		 * @param[in] end iterator to the end of the upper bounds sequence.
		 */
		template <class Iterator>
		void set_ub(Iterator start, Iterator end)
		{
			if (std::distance(start,end) != std::distance(m_lb.begin(),m_lb.end())) {
				pagmo_throw(value_error,"invalid bounds dimension in set_ub()");
			}
			verify_bounds(m_lb.begin(),m_lb.end(),start,end);
			std::copy(start,end,m_ub.begin());
			// Normalise bounds.
			normalise_bounds();
		}
		/// Upper bounds setter from raw array.
		/**
		 * Will fail if N is different from global problem dimension or if at least one upper bound is less than the
		 * corresponding lower bound.
		 *
		 * @param[in] v upper bounds array.
		 */
		template <std::size_t N>
		void set_ub(const double (&v)[N])
		{
			if (N != m_lb.size()) {
				pagmo_throw(value_error,"invalid bounds dimension in set_ub()");
			}
			verify_bounds(m_lb.begin(),m_lb.end(),v,v + N);
			std::copy(v,v + N,m_ub.begin());
			// Normalise bounds.
			normalise_bounds();
		}
		//@}
		/** @name Properties getters.*/
		//@{
		unsigned int get_fevals() const;
		unsigned int get_cevals() const;
		size_type get_dimension() const;
		size_type get_i_dimension() const;
		f_size_type get_f_dimension() const;
		c_size_type get_c_dimension() const;
		c_size_type get_ic_dimension() const;
		const std::vector<double>& get_c_tol() const;
		double get_diameter() const;
		virtual std::string get_name() const;
		//@}
		constraint_vector compute_constraints(const decision_vector &) const;
		void compute_constraints(constraint_vector &, const decision_vector &) const;
		bool compare_constraints(const constraint_vector &, const constraint_vector &) const;
		bool test_constraint(const constraint_vector &, const c_size_type &) const;
		bool feasibility_x(const decision_vector &) const;
		bool feasibility_c(const constraint_vector &) const;
		/// Clone method.
		/**
		 * Provided that the derived problem implements properly the copy constructor, virtually all implementations of this method will
		 * look like this:
@verbatim
return base_ptr(new derived_problem(*this));
@endverbatim
		 *
		 * @return problem::base_ptr to a copy of this.
		 */
		virtual base_ptr clone() const = 0;
		std::string human_readable() const;
		virtual std::string human_readable_extra() const;
		bool operator==(const base &) const;
		bool operator!=(const base &) const;
		bool is_compatible(const base &) const;
		bool compare_x(const decision_vector &, const decision_vector &) const;
		bool verify_x(const decision_vector &) const;
		bool compare_fc(const fitness_vector &, const constraint_vector &, const fitness_vector &, const constraint_vector &) const;
		virtual void pre_evolution(population &) const;
		virtual void post_evolution(population &) const;
	protected:
		virtual bool equality_operator_extra(const base &) const;
		virtual void compute_constraints_impl(constraint_vector &, const decision_vector &) const;
		virtual bool compare_constraints_impl(const constraint_vector &, const constraint_vector &) const;
		virtual bool compare_fc_impl(const fitness_vector &, const constraint_vector &, const fitness_vector &, const constraint_vector &) const;
		void estimate_sparsity(const decision_vector &, int& lenG, std::vector<int>& iGfun, std::vector<int>& jGvar) const;
		void estimate_sparsity(int& lenG, std::vector<int>& iGfun, std::vector<int>& jGvar) const;
	public:
		virtual void set_sparsity(int& lenG, std::vector<int>& iGfun, std::vector<int>& jGvar) const;
		/** @name Objective function and fitness handling.
		 * Methods used to calculate and compare fitnesses.
		 */
		//@{
		fitness_vector objfun(const decision_vector &) const;
		void objfun(fitness_vector &, const decision_vector &) const;
		bool compare_fitness(const fitness_vector &, const fitness_vector &) const;
		void reset_caches() const;
	public:
		const std::vector<constraint_vector>& get_best_c(void) const;
		const std::vector<decision_vector>& get_best_x(void) const;
		const std::vector<fitness_vector>& get_best_f(void) const;
		void set_best_x(const std::vector<decision_vector>&);
	protected:
		virtual bool compare_fitness_impl(const fitness_vector &, const fitness_vector &) const;
		/// Objective function implementation.
		/**
		 * Takes a pagmo::decision_vector x as input and writes its pagmo::fitness_vector to f. This function is not to be called directly,
		 * it is invoked by objfun() after a series of safety checks is performed on x and f.
		 *
		 * @param[out] f fitness vector into which x's fitness will be written.
		 * @param[in] x decision vector whose fitness will be calculated.
		 */
		virtual void objfun_impl(fitness_vector &f, const decision_vector &x) const = 0;
		//@}
	private:
		void normalise_bounds();
		// Construct from iterators.
		template <class Iterator1, class Iterator2>
		void construct_from_iterators(Iterator1 start1, Iterator1 end1, Iterator2 start2, Iterator2 end2)
		{
			m_lb.insert(m_lb.end(),start1,end1);
			m_ub.insert(m_ub.end(),start2,end2);
			if (m_lb.size() != m_ub.size() || m_lb.size() == 0) {
				pagmo_throw(value_error,"null or inconsistent dimension(s) for upper/lower bounds while constructing problem");
			}
			verify_bounds(m_lb.begin(),m_lb.end(),m_ub.begin(),m_ub.end());
		}
		// Verify upper/lower bounds. This must be called only after having made sure that the iterator distances
		// are consistent.
		template <class Iterator1, class Iterator2>
		static void verify_bounds(Iterator1 start1, Iterator1 end1, Iterator2 start2, Iterator2 end2)
		{
			for (; start1 != end1 && start2 != end2; ++start1, ++start2) {
				if (*start1 > *start2) {
					pagmo_throw(value_error,"lower bound is greater than upper bound");
				}
			}
		}
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & const_cast<size_type &>(m_i_dimension);
			ar & const_cast<f_size_type &>(m_f_dimension);
			ar & const_cast<c_size_type &>(m_c_dimension);
			ar & const_cast<c_size_type &>(m_ic_dimension);
			ar & m_lb;
			ar & m_ub;
			ar & const_cast<std::vector<double> &>(m_c_tol);
			ar & m_decision_vector_cache_f;
			ar & m_fitness_vector_cache;
			ar & m_decision_vector_cache_c;
			ar & m_constraint_vector_cache;
			ar & m_tmp_f1;
			ar & m_tmp_f2;
			ar & m_tmp_c1;
			ar & m_tmp_c2;
			ar & m_best_x;
			ar & m_best_f;
			ar & m_best_c;
			ar & m_fevals;
			ar & m_cevals;
		}

		// Data members.
		// Size of the integer part of the problem.
		const size_type				m_i_dimension;
		// Size of the fitness vector.
		const f_size_type			m_f_dimension;
		// Global constraints dimension.
		const c_size_type			m_c_dimension;
		// Inequality constraints dimension
		const c_size_type			m_ic_dimension;
		// Lower bounds.
		decision_vector				m_lb;
		// Upper bounds.
		decision_vector				m_ub;
		// Tolerance for constraints analysis.
		const std::vector<double>   m_c_tol;
		// Decision vector cache for fitness.
		mutable decision_vector_cache_type	m_decision_vector_cache_f;
		// Fitness vector cache.
		mutable fitness_vector_cache_type	m_fitness_vector_cache;
		// Decision vector cache for constraints.
		mutable decision_vector_cache_type	m_decision_vector_cache_c;
		// Constraint vector cache.
		mutable constraint_vector_cache_type	m_constraint_vector_cache;
		// Temporary storage used during decision_vector comparisons.
		mutable fitness_vector			m_tmp_f1;
		mutable fitness_vector			m_tmp_f2;
		// Temporary storage used during constraints satisfaction testing and constraints comparison.
		mutable constraint_vector		m_tmp_c1;
		mutable constraint_vector		m_tmp_c2;

		// Best known vectors
		std::vector<decision_vector> m_best_x;
		std::vector<fitness_vector> m_best_f;
		std::vector<constraint_vector> m_best_c;

		// Number of function and constraints evaluations
		mutable unsigned int                    m_fevals;
		mutable unsigned int                    m_cevals;
};

std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const base &);

}
}

BOOST_SERIALIZATION_ASSUME_ABSTRACT(pagmo::problem::base)

#endif
