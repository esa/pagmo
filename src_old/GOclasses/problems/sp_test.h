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

// 09/09/09 Created by Dario Izzo.

#ifndef PAGMO_PROBLEM_SP_TEST_H
#define PAGMO_PROBLEM_SP_TEST_H

#include <string>
#include <vector>

#include "../../config.h"
#include "base.h"

namespace pagmo
{
namespace problem {

/// Stochastic Programming Test Problem: Inventory Model
/**
 * This simple test problem represent an inventory problem and is discussed e.g. in A Tutorial on Stochastic Programming
 * by Shapiro.
 */
class __PAGMO_VISIBLE inventory: public base
{
	public:
		inventory(int sample_size);
		virtual inventory *clone() const {
			return new inventory(*this);
		}
		virtual void pre_evolution(population &) const;
		virtual std::string id_object() const {
			return id_name();
		}
	private:
		virtual double objfun_(const std::vector<double> &) const;
		mutable std::vector<double> d;
		size_t m_sample_size;
};

}
}

#endif
