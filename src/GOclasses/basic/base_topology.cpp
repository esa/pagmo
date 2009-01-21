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

#include "archipelago.h"
#include "base_topology.h"
#include "island.h"
#include "population.h"

base_topology::base_topology() {}

base_topology::~base_topology() {}

Population &base_topology::get_pop(island &isl)
{
	return isl.m_pop;
}

const Population &base_topology::get_pop(const island &isl)
{
	return isl.m_pop;
}

archipelago *base_topology::get_arch(island &isl)
{
	return isl.m_a;
}

const archipelago *base_topology::get_arch(const island &isl)
{
	return isl.m_a;
}
