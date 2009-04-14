/*****************************************************************************
 *   Copyright (C) 2008, 2009 Advanced Concepts Team (European Space Agency) *
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

// 09/01/2009: Initial version by Francesco Biscani.

#ifndef PAGMO_IHS_ALGORITHM_H
#define PAGMO_IHS_ALGORITHM_H

#include <iostream>

#include "../../config.h"
#include "../basic/population.h"
#include "go_algorithm.h"

/// Improved harmony search algorithm.
/**
 * See http://en.wikipedia.org/wiki/Harmony_search
 */
class __PAGMO_VISIBLE ihs_algorithm: public go_algorithm {
	public:
		ihs_algorithm(int, const double &, const double &, const double &, const double &, const double &);
		virtual Population evolve(const Population &) const;
		virtual ihs_algorithm *clone() const {return new ihs_algorithm(*this);}
		virtual std::string id_object() const {return id_name(); }
	private:
		virtual void log(std::ostream &) const;
		// Number of generations.
		const size_t	m_gen;
		// Rate of choosing from memory (i.e., population).
		const double	m_phmcr;
		// Minimum pitch adjustment rate.
		const double	m_ppar_min;
		// Maximum pitch adjustment rate.
		const double	m_ppar_max;
		// Mininum distance bandwidth.
		const double	m_bw_min;
		// Maximum distance bandwidth.
		const double	m_bw_max;
};

#endif
