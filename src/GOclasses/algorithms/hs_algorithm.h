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

#ifndef PAGMO_HS_ALGORITHM_H
#define PAGMO_HS_ALGORITHM_H

#include "../../config.h"
#include "../basic/population.h"
#include "go_algorithm.h"

/// Harmony search algorithm.
/**
 * See http://en.wikipedia.org/wiki/Harmony_search.
 */
class __PAGMO_VISIBLE hs_algorithm: public go_algorithm {
	public:
		hs_algorithm(int, const double &, const double &, const double &);
		virtual Population evolve(const Population &) const;
		virtual hs_algorithm *clone() const {return new hs_algorithm(*this);}
	private:
		// Number of generations.
		const size_t	m_gen;
		// Rate of choosing from memory (i.e., population).
		const double	m_phmcr;
		// Pitch adjustment rate.
		const double	m_ppar;
		// 'Distance bandwidth', i.e. max change for pitch adjustment.
		const double	m_bw;
};

#endif
