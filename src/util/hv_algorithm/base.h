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

#ifndef PAGMO_UTIL_HV_ALGORITHM_BASE_H
#define PAGMO_UTIL_HV_ALGORITHM_BASE_H

#include <iostream>
#include <string>
#include <typeinfo>
#include <boost/lexical_cast.hpp>

#include "../../config.h"
#include "../../exceptions.h"
#include "../../serialization.h"
#include "../../types.h"

namespace pagmo { namespace util {
/// Hypervolume algorithm namespace.
/**
 * This namespace contains all the algorithms implemented for the purpose of calculating the hypervolume indicator.
 */
namespace hv_algorithm {

/// Base hypervolume algorithm class.
class base;

/// Base hypervolume algorithm class
typedef boost::shared_ptr<base> base_ptr;

/// Base hypervolume class.
/**
 * This class represents the abstract hypervolume algorithm used for computing
 * the hypervolume indicator (also known as Lebesgue measure, or S-metric).
 *
 * Every hypervolume algorithm that extends this base class implement the following methods:
 * - compute() which takes a list of points and a reference point.
 * - verify_before_compute() which raises an exception for special cases of misuse of this method (e.g. method working only for some number of dimensions).
 *
 * @author Krzysztof Nowak (kn@kiryx.net)
 */
class __PAGMO_VISIBLE base
{
	public:

		/// Compute method
		/**
		 * This method computes the hypervolume.
		 * It accepts a list of points as an input, and the distinguished "reference point".
		 * Hypervolume is then computed as a joint hypervolume of hypercubes, generated pairwise with the reference point and each point from the set.
		 *
		 * @param[in] points - vector of fitness_vectors for which the hypervolume is computed
		 * @param[in] r_point - distinguished "reference point".
		 */
		virtual double compute(std::vector<fitness_vector> &points, const fitness_vector &r_point) = 0;

		virtual double exclusive(const unsigned int p_idx, std::vector<fitness_vector> &points, const fitness_vector &r_point);

		virtual unsigned int least_contributor(std::vector<fitness_vector> &points, const fitness_vector &r_point);

		virtual unsigned int greatest_contributor(std::vector<fitness_vector> &points, const fitness_vector &r_point);

		/// Verification of input
		/**
		 * This method serves as a verification method.
		 * Not every algorithm is suited of every type of problem.
		 *
		 * @param[in] points - vector of fitness_vectors for which the hypervolume is computed
		 * @param[in] r_point - distinguished "reference point".
		 */
		virtual void verify_before_compute(const std::vector<fitness_vector> &points, const fitness_vector &r_point) = 0;

		/// Clone method.
		/**
		 * @return pagmo::util::hv_algorithm::base_ptr to a copy of this.
		 */
		virtual base_ptr clone() const = 0;

		virtual std::string get_name() const;
		virtual ~base();

	protected:

		void assert_maximal_reference_point(const std::vector<fitness_vector> &points, const fitness_vector &r_point);

		virtual unsigned int extreme_contributor(std::vector<fitness_vector> &points, const fitness_vector &r_point, bool (*)(double, double));

		// comparison functions for the least and the greatest contributor methods
		static bool cmp_least(const double, const double);
		static bool cmp_greatest(const double, const double);

		/// Compute volume between two points
		/**
		 * Calculates the volume between points a and b (as defined for n-dimensional Euclidean spaces).
		 *
		 * @param[in] a first point defining the hypercube
		 * @param[in] b second point defining the hypercube
		 * @param[in] dim_bound dimension boundary for the volume. If equal to 0, then compute the volume of whole vector. Any positive number limits the computation from dimension 0 to dim_bound INCLUSIVE.
		 *
		 * @return volume of hypercube defined by points a and b
		 */
		inline double volume_between(const fitness_vector &a, const fitness_vector &b, unsigned int dim_bound = 0) const {
			if (dim_bound == 0) { 
				dim_bound = a.size();
			}
			double volume = 1.0;
			for (fitness_vector::size_type idx = 0; idx < dim_bound ; ++idx) {
				volume *= (a[idx] - b[idx]);
			}
			return (volume < 0 ? -volume : volume);
		}

		/// Compute volume between two points
		/**
		 * Calculates the volume between points a and b (as defined for n-dimensional Euclidean spaces).
		 *
		 * @param[in] a first point defining the hypercube
		 * @param[in] b second point defining the hypercube
		 * @param[in] size dimension of the vectors.
		 *
		 * @return volume of hypercube defined by points a and b
		 */
		inline double volume_between(double* a, double* b, unsigned int size) const {
			double volume = 1.0;
			while(size--) {
				volume *= (b[size] - a[size]);

			}
			return (volume < 0 ? -volume : volume);
		}

		// Domination results of the 'dom_cmp' methods
		enum {
			DOM_CMP_B_DOMINATES_A = 1, // second argument dominates the first one
			DOM_CMP_A_DOMINATES_B = 2, // first argument dominates the second one
			DOM_CMP_A_B_EQUAL = 3, // both points are equal
			DOM_CMP_INCOMPARABLE = 4 // points are incomparable
		};

		/// dominance comparison method
		/**
		 * Establishes the domination relationship between two points.
		 *
		 * returns DOM_CMP_B_DOMINATES_A if point 'b' DOMINATES point 'a'
		 * returns DOM_CMP_A_DOMINATES_B if point 'a' DOMINATES point 'b'
		 * returns DOM_CMP_A_B_EQUAL if point 'a' IS EQUAL TO 'b'
		 * returns DOM_CMP_INCOMPARABLE otherwise
		 */
		inline int dom_cmp(const fitness_vector &a, const fitness_vector &b, unsigned int dim_bound = 0) const {
			if (dim_bound == 0) {
				dim_bound = a.size();
			}
			for(fitness_vector::size_type i = 0; i < dim_bound ; ++i) {
				if (a[i] > b[i]) {
					for(fitness_vector::size_type j = i + 1; j < dim_bound ; ++j) {
						if (a[j] < b[j]) {
							return DOM_CMP_INCOMPARABLE;
						}
					}
					return DOM_CMP_B_DOMINATES_A;
				}
				else if (a[i] < b[i]) {
					for(fitness_vector::size_type j = i + 1 ; j < dim_bound ; ++j) {
						if (a[j] > b[j]) {
							return DOM_CMP_INCOMPARABLE;
						}
					}
					return DOM_CMP_A_DOMINATES_B;
				}
			}
			return DOM_CMP_A_B_EQUAL;
		}

		/// dominance comparison method
		/**
		 * Establishes the domination relationship between two points (overloaded for double*);
		 *
		 * returns DOM_CMP_B_DOMINATES_A if point 'b' DOMINATES point 'a'
		 * returns DOM_CMP_A_DOMINATES_B if point 'a' DOMINATES point 'b'
		 * returns DOM_CMP_A_B_EQUAL if point 'a' IS EQUAL TO 'b'
		 * returns DOM_CMP_INCOMPARABLE otherwise
		 */
		inline int dom_cmp(double *a, double* b, unsigned int size) const {
			for(fitness_vector::size_type i = 0; i < size ; ++i) {
				if (a[i] > b[i]) {
					for(fitness_vector::size_type j = i + 1; j < size ; ++j) {
						if (a[j] < b[j]) {
							return DOM_CMP_INCOMPARABLE;
						}
					}
					return DOM_CMP_B_DOMINATES_A;
				}
				else if (a[i] < b[i]) {
					for(fitness_vector::size_type j = i + 1 ; j < size ; ++j) {
						if (a[j] > b[j]) {
							return DOM_CMP_INCOMPARABLE;
						}
					}
					return DOM_CMP_A_DOMINATES_B;
				}
			}
			return DOM_CMP_A_B_EQUAL;
		}

	private:

		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int) {
			(void)ar;
		}
};

///Fitness vector comparator class
/**
 * This is a helper class that allows for the generation of comparator objects.
 * Many hypervolume algorithms use comparator functions for sorting, or data structures handling.
 * In most cases the difference between the comparator functions differ either by the dimension number, or the inequality sign ('>' or '<').
 * We provide a general comparator class for that purpose.
 */
class fitness_vector_cmp {

	public:
		fitness_vector_cmp(int dim, char cmp_type);
	
		///Overloaded operator()
		/**
		 * Overloading call operator is required for all sorting and data structure key comparators in stl library.
		 *
		 * @param[in] lhs fitness_vector on the left hand side
		 * @param[in] rhs fitness_vector on the right hand side
		 *
		 * @return Boolean variable stating whether given expression is true for fitness_vectors.
		 */
		inline bool operator()(const fitness_vector &lhs, const fitness_vector &rhs) {
			return (*m_cmp_obj)(lhs,rhs);
		}
	private:
		struct cmp_fun {
			int m_dim;
			cmp_fun(int dim) : m_dim(dim) { }
			virtual ~cmp_fun() { };
			/// virtual operator() - It is never called anyway, so we could have gone with pure virtual, yet then we would not be able to use inline.
			virtual inline bool operator()(const fitness_vector &lhs, const fitness_vector &rhs) {
				return lhs[0] < rhs[0];
			}
		};

		struct cmp_le : cmp_fun {	
			cmp_le(int dim) : cmp_fun(dim) { }
			inline bool operator()(const fitness_vector &lhs, const fitness_vector &rhs) {
				return lhs[m_dim] < rhs[m_dim];
			}
		};

		struct cmp_ge : cmp_fun {	
			cmp_ge(int dim) : cmp_fun(dim) { }
			inline bool operator()(const fitness_vector &lhs, const fitness_vector &rhs) {
				return lhs[m_dim] > rhs[m_dim];
			}
		};

		struct boost::shared_ptr<cmp_fun> m_cmp_obj;
};

} } }

BOOST_SERIALIZATION_ASSUME_ABSTRACT(pagmo::util::hv_algorithm::base);

#endif
