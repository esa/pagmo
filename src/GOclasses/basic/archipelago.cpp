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

#include "archipelago.h"

archipelago::archipelago(const archipelago &a):m_container(a.m_container) {}

size_t archipelago::size() const
{
	return m_container.size();
}

void archipelago::push_back(const value_type &i)
{
	m_container.push_back(i);
	m_container.back().set_archipelago(this);
}

// const island &archipelago::operator[](const size_t &n) const
// {
// 	return m_container[n];
// }

// island &archipelago::operator[](const size_t &n)
// {
// 	return m_container[n];
// }


