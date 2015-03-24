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

#ifndef PAGMO_MIGRATION_BASE_POLICY_H
#define PAGMO_MIGRATION_BASE_POLICY_H

#include <iostream>
#include <string>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"

namespace pagmo {

/// Migration policies namespace.
/**
 * This namespace contains selection/replacement policies used during migration in the archipelago class.
 */
namespace migration {

/// Type of migration rate.
/**
 * Used both by selection and replacement migration policies.
 */
enum rate_type
{
	/// Migration rate is interpreted as the absolute number of individuals to migrate.
	absolute = 0,
	/// Migration rate is interpreted as the fraction of individuals to migrate with respect to the orign/destination population.
	fractional = 1
};

/// Base migration class.
/**
 * Embeds two properties used both in selection and replacement policies:
 * - a migration rate type, which can be either absolute or fractional;
 * - a migration rate, which represents:
 *   - the absolute number of individuals selected/replaced from/in the interested population (absolute migration rate type),
 *   - the fraction of individuals selected/replaced from/in the interested population (fractional migration rate type).
 *
 * The get_n_individuals() method returns the number of individuals to be selected/replaced in the given population, according to the
 * properties above.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 * @author Marek Ruciński (marek.rucinski@gmail.com)
 */
class __PAGMO_VISIBLE base
{
	public:
		base(const double &, rate_type);
		population::size_type get_n_individuals(const population &) const;
		std::string human_readable() const;
		virtual ~base();
	protected:
		virtual std::string human_readable_extra() const;
	protected:
		/// Iota function, usefull to fill iterator range with increasing values.
		/**
		 * This function will fill the contents in the range [first,last) with increasing values, starting from *first = value.
		 *
		 * @param[in] first start of the iterator range.
		 * @param[in] last end of the iterator range.
		 * @param[in] value initial value.
		 */
		template <class ForwardIterator, class T>
		static void iota(ForwardIterator first, ForwardIterator last, T value)
		{
			for (; first != last; ++first, ++value)
				*first = value;
		}
		/// Migration rate.
		/**
		 * It will be interpreted as an integer in case of absolute rate migration type, as a floating-point value
		 * in case of fractional migration type.
		 */
		double		m_rate;
		/// Migration rate type.
		rate_type	m_type;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
				ar & m_rate;
				ar & m_type;
		}
};

std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const base &);

} }

#endif
