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

#ifndef PAGMO_ALGORITHM_BASE_H
#define PAGMO_ALGORITHM_BASE_H

#include <iostream>
#include <string>
#include <typeinfo>
#include <boost/shared_ptr.hpp>

#include "../config.h"
#include "../population.h"
#include "../rng.h"
#include "../serialization.h"

namespace pagmo
{
/// Algorithm namespace.
/**
 * This namespace contains all the algorithms implemented in PaGMO.
 */
namespace algorithm {

/// Base algorithm class.
class base;

/// Alias for shared pointer to base algorithm.
typedef boost::shared_ptr<base> base_ptr;

/// Base algorithm class.
/**
 * All algorithms implemented in PaGMO must derive from this base class. This base class provides each algorithm with one pagmo::rng_double
 * and one pagmo::rng_uint32 random number generators. Each algorithm must implement the base::evolve() method.
 *
 * \section Serialization
 * The algorithm classes are serialized for the purpose of being transmitted over a distributed environment.
 * Serializing a derived algorithm requires that the needed serialization libraries be declared in their header of the derived class. 
 * Virtually all the derived algorithm classes need to have the following declared in the header files:
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
 * In order to be able to identify the dervied class from a base_pointer the derived class needs to be registered. This is done by registering the class in the "pagmo/src/algorithms.h", in the REGISTER_ALGORITHM_SERIALIZATIONS() routine.
 * Notes: 
 * - "const" attributes need to be cast as constants in the serialize method using const_cast
 * - attributes that that are not primitives, need be a serialized type as well
 * - pointers to primitives cannot be serialized (in this case one can split the serialize method into save/load methods and store the values, that the pointers refer to, into temporary variables which are serialized insted - see boost serialize documentation on the topic if needed)
 *
 * @author Francesco Biscani (bluescarni@gmail.com) 
 */
class __PAGMO_VISIBLE base
{    
	public:
		base();
		/// Evolve method.
		/**
		 * The purpose of this method is to take a pagmo::population as input and evolve it towards the solution of the problem.
		 *
		 * @param[in,out] p population to be evolved.
		 */
		virtual void evolve(population &p) const = 0;
		/// Clone method.
		/**
		 * Provided that the derived algorithm implements properly the copy constructor, virtually all implementations of this method will
		 * look like this:
@verbatim
return base_ptr(new derived_algorithm(*this));
@endverbatim
		 *
		 * @return algorithm::base_ptr to a copy of this.
		 */
		virtual base_ptr clone() const = 0;
		virtual ~base();
		std::string human_readable() const;
		virtual std::string get_name() const;
		virtual std::string human_readable_extra() const;
		/// Setter-Getter for protected m_screen_output data
		void set_screen_output(const bool p);
		bool get_screen_output() const;

		/// Resets the seed of the internal rngs using a user-provided seed
		void reset_rngs(const unsigned int) const;

	protected:
		/// Indicates to the derived class whether to print stuff on screen
		bool m_screen_output;
		/// Random number generator for double-precision floating point values.
		mutable rng_double	m_drng;
		/// Random number generator for unsigned integer values.
		mutable rng_uint32	m_urng;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & m_drng;
			ar & m_urng; 
			ar & m_screen_output;
		}
	protected:
		/// A counter for the number of function evaluations
		mutable unsigned int m_fevals;
};

std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const base &);

}
}

BOOST_SERIALIZATION_ASSUME_ABSTRACT(pagmo::algorithm::base)

#endif
