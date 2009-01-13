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

// 12/01/2009: Initial version by Francesco Biscani.

#ifndef PAGMO_BASE_TOPOLOGY_H
#define PAGMO_BASE_TOPOLOGY_H

#include "../../../config.h"
#include "population.h"

class archipelago;

class __PAGMO_VISIBLE base_topology {
		friend class archipelago;
	public:
		base_topology();
		virtual ~base_topology();
		virtual void pre_evolution(const Population &) = 0;
		virtual void post_evolution(const Population &) = 0;
		virtual base_topology *clone() const = 0;
		archipelago &arch();
	protected:
		void set_archipelago(archipelago *);
		archipelago	*m_arch;
};

#endif
