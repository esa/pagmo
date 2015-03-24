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
 * the hypervolume indicator (also known as Lebesgue measure, or S-metric), and other
 * measures that can derive from it, e.g. exclusive contribution by given point
 *
 * The are few public methods available:
 *
 * - 'compute' - returns the total hypervolume of the set of points
 * - 'exclusive' - returns the exclusive volume contributed by given point
 * - 'least_contributor' - returns an index of the point contributing the least volume
 * - 'greatest_contributor' - returns an index of the point contributing the most volume
 * - 'contributions' - returns a vector of exclusive contributions by each of the points.
 *
 * Additionally a private method 'base::extreme_contributor' can be overloaded:
 * - 'extreme_contributor' - returns an index of a single individual that contributes either the least or the greatest
 *  amount of the volume. The distinction between the extreme contributors is made using a comparator function.
 *  Purpose of this method is to avoid repeating a similar code for the least and the greatest contributor methods.
 *  In many cases it's just a matter of a single alteration in a comparison sign '<' to '>'. See concrete example here for more details.
 *
 * The following base class provides an interface for any hv_algorithm that may be added.
 * The crucial method to implement is the 'compute' method, as the remaining methods can be derived from it.
 *
 * Base class assumes that the hv_algorithm implements the 'compute' method, and employs a naive approach to provide other functionalities:
 *
 * 'base::exclusive' method relies on 'compute' method, by computing the hypervolume twice (e.g. ExclusiveHV = HV(P) - HV(P/{p}))
 * 'base::contributions' method relies on 'compute' method as well, by computing the exclusive volume for each point using approach above.
 * 'base::extreme_contributor' (private method) relies on the 'base::contributions' method in order to elicit the correct extreme individual.
 * 'base::least_contributor' and 'base::greatest_contributor' methods rely on 'base::extreme_contributor' method by providing the correct comparator.
 *
 * Thanks to that, any newly implemented hypervolume algorithm that overloads the 'compute' method, gets the functionalities above as well.
 * It is often the case that the algorithm may provide a better solution for each of the features above, e.g. overloading the 'base::extreme_contributor'
 * method with an efficient solution, will automatically speed up the 'least_contributor' and the 'greatest_contributor' methods as well.
 *
 * Additionally, any newly implemented hypervolume algorithm should overload the 'base::verify_before_compute' method in order to prevent
 * the computation for the incompatible data.
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
	virtual double compute(std::vector<fitness_vector> &points, const fitness_vector &r_point) const = 0;
	virtual double exclusive(const unsigned int, std::vector<fitness_vector> &, const fitness_vector &) const;
	virtual unsigned int least_contributor(std::vector<fitness_vector> &, const fitness_vector &) const;
	virtual unsigned int greatest_contributor(std::vector<fitness_vector> &, const fitness_vector &) const;
	virtual std::vector<double> contributions(std::vector<fitness_vector> &, const fitness_vector &) const;

	/// Verification of input
	/**
	 * This method serves as a verification method.
	 * Not every algorithm is suited of every type of problem.
	 *
	 * @param[in] points - vector of fitness_vectors for which the hypervolume is computed
	 * @param[in] r_point - distinguished "reference point".
	 */
	virtual void verify_before_compute(const std::vector<fitness_vector> &points, const fitness_vector &r_point) const = 0;

	/// Clone method.
	/**
	 * @return pagmo::util::hv_algorithm::base_ptr to a copy of this.
	 */
	virtual base_ptr clone() const = 0;

	virtual std::string get_name() const;
	virtual ~base();

protected:
	void assert_minimisation(const std::vector<fitness_vector> &, const fitness_vector &) const;

	virtual unsigned int extreme_contributor(std::vector<fitness_vector> &, const fitness_vector &, bool (*)(double, double)) const;

	// comparison functions for the least and the greatest contributor methods
	static bool cmp_least(const double, const double);
	static bool cmp_greatest(const double, const double);

public:
	static double volume_between(const fitness_vector &, const fitness_vector &, unsigned int = 0);

protected:
	// Domination results of the 'dom_cmp' methods
	enum {
		DOM_CMP_B_DOMINATES_A = 1, // second argument dominates the first one
		DOM_CMP_A_DOMINATES_B = 2, // first argument dominates the second one
		DOM_CMP_A_B_EQUAL = 3, // both points are equal
		DOM_CMP_INCOMPARABLE = 4 // points are incomparable
	};

	static double volume_between(double*, double*, unsigned int);
	static int dom_cmp(double*, double*, unsigned int);
	static int dom_cmp(const fitness_vector &, const fitness_vector &, unsigned int = 0);

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int) {
		(void)ar;
	}
};

/// Fitness vector comparator class
/**
 * This is a helper class that allows for the generation of comparator objects.
 * Many hypervolume algorithms use comparator functions for sorting, or data structures handling.
 * In most cases the difference between the comparator functions differ either by the dimension number, or the inequality sign ('>' or '<').
 * We provide a general comparator class for that purpose.
 */
class fitness_vector_cmp
{
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
	inline bool operator()(const fitness_vector &lhs, const fitness_vector &rhs)
	{
		return (*m_cmp_obj)(lhs,rhs);
	}
private:
	struct cmp_fun
	{
		int m_dim;
		cmp_fun(int dim) : m_dim(dim) { }
		virtual ~cmp_fun() { };
		/// virtual operator() - It is never called anyway, so we could have gone with pure virtual, yet then we would not be able to use inline.
		virtual inline bool operator()(const fitness_vector &lhs, const fitness_vector &rhs)
		{
			return lhs[0] < rhs[0];
		}
	};

	struct cmp_le : cmp_fun
	{
		cmp_le(int dim) : cmp_fun(dim) { }
		inline bool operator()(const fitness_vector &lhs, const fitness_vector &rhs)
		{
			return lhs[m_dim] < rhs[m_dim];
		}
	};

	struct cmp_ge : cmp_fun
	{
		cmp_ge(int dim) : cmp_fun(dim) { }
		inline bool operator()(const fitness_vector &lhs, const fitness_vector &rhs)
		{
			return lhs[m_dim] > rhs[m_dim];
		}
	};

	boost::shared_ptr<cmp_fun> m_cmp_obj;
};

} } }

BOOST_SERIALIZATION_ASSUME_ABSTRACT(pagmo::util::hv_algorithm::base)

#endif
