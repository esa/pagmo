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

// 04/01/2009: Initial version by Francesco Biscani.

#include <typeinfo>

#include "../../exceptions.h"
#include "../problems/GOproblem.h"
#include "../algorithms/go_algorithm.h"
#include "archipelago.h"
#include "base_topology.h"
#include "no_topology.h"
#include "island.h"

archipelago::archipelago(const GOProblem &p):m_gop(p.clone()),m_top(new no_topology()) {}

archipelago::archipelago(const GOProblem &p, const base_topology &t):m_gop(p.clone()),m_top(t.clone()) {}

archipelago::archipelago(const GOProblem &p, const go_algorithm &a, int N, int M):m_gop(p.clone()),m_top(new no_topology())
{
	for (int i = 0; i < N; ++i) {
		push_back(island(p,a,M));
	}
}

archipelago::archipelago(const GOProblem &p, const base_topology &t, const go_algorithm &a, int N, int M):m_gop(p.clone()),m_top(t.clone())
{
	for (int i = 0; i < N; ++i) {
		push_back(island(p,a,M));
	}
}

archipelago::archipelago(const archipelago &a):m_container(a.m_container),m_gop(a.m_gop->clone()),m_top(a.m_top->clone())
{
	const iterator it_f = end();
	for (iterator it = begin(); it != it_f; ++it) {
		it->set_archipelago(this);
	}
}

archipelago &archipelago::operator=(const archipelago &a)
{
	// We want to guard against assignment here because the topology might access the info on the island
	// list while the list itself is being modified (please note that if it weren't for the topology this would be fine,
	// since the single island take care of themselves against concurrent read/write access through internal
	// mutexes).
	join();
	if (this != &a) {
		if (typeid(*m_gop) != typeid(*a.m_gop)) {
			pagmo_throw(type_error, "problem types are not compatible while assigning archipelago");
		}
		m_container = a.m_container;
		const iterator it_f = end();
		for (iterator it = begin(); it != it_f; ++it) {
			it->set_archipelago(this);
		}
		m_gop.reset(a.m_gop->clone());
		m_top.reset(a.m_top->clone());
	}
	return *this;
}

const base_topology &archipelago::topology() const
{
	return *m_top;
}

void archipelago::set_topology(const base_topology &t)
{
	join();
	m_top.reset(t.clone());
}

island &archipelago::operator[](int n)
{
	// NOTE: maybe this is not really needed, island should be able to take care of itself
	// without external protection.
	join();
	return *it_from_index<iterator>(n);
}

size_t archipelago::size() const
{
	return m_container.size();
}

void archipelago::check_island(const island &isl) const
{
	if (typeid(isl.problem()) != typeid(problem())) {
		pagmo_throw(type_error, "island's problem type is not compatible with archipelago's problem type");
	}
}

size_t archipelago::island_index(const island *isl) const
{
	const const_iterator it_f = end();
	size_t retval = 0;
	for (const_iterator it = begin(); it != it_f; ++it, ++retval) {
		if (&(*it) == isl) {
			return retval;
		}
	}
	pagmo_throw(index_error,"could not find island index");
}

void archipelago::push_back(const island &isl)
{
	join();
	check_island(isl);
	m_container.push_back(isl);
	m_container.back().set_archipelago(this);
}

void archipelago::insert(int n, const island &isl)
{
	join();
	check_island(isl);
	iterator it;
	// Emulate Python lists' behaviour of inserting at the beginning if the negative index
	// is out of range, at the end if the positive index is out of range.
	try {
		it = it_from_index<iterator>(n);
	} catch (const index_error &) {
		if (n >= 0) {
			it = m_container.begin();
		} else {
			it = m_container.end();
		}
	}
	m_container.insert(it,isl);
	--it;
	it->set_archipelago(this);
}

const GOProblem &archipelago::problem() const
{
	return *m_gop;
}

void archipelago::join() const
{
	const const_iterator it_f = m_container.end();
	for (const_iterator it = m_container.begin(); it != it_f; ++it) {
		it->join();
	}
}

bool archipelago::busy() const
{
	const const_iterator it_f = m_container.end();
	for (const_iterator it = m_container.begin(); it != it_f; ++it) {
		if (it->busy()) {
			return true;
		}
	}
	return false;
}

void archipelago::evolve(int n)
{
	if (busy()) {
		pagmo_throw(runtime_error,"cannot start evolution while evolving");
	}
	const iterator it_f = m_container.end();
	for (iterator it = m_container.begin(); it != it_f; ++it) {
		it->evolve(n);
	}
}

void archipelago::evolve_t(const size_t &t)
{
	if (busy()) {
		pagmo_throw(runtime_error,"cannot start evolution while evolving");
	}
	const iterator it_f = m_container.end();
	for (iterator it = m_container.begin(); it != it_f; ++it) {
		it->evolve_t(t);
	}
}

std::ostream &operator<<(std::ostream &s, const archipelago &a) {
	s << "Problem type:    " << a.m_gop->id_name() << "\n\n";
	const archipelago::const_iterator it_f = a.m_container.end();
	size_t i = 0;
	for (archipelago::const_iterator it = a.m_container.begin(); it != it_f; ++it) {
		s << "Island #:        " << i << '\n';
		s << "ID:              " << it->id() << '\n';
		s << "Population size: " << it->size() << '\n';
		s << "Evolution time:  " << it->evo_time() << '\n';
		s << "Algorithm type:  " << it->algorithm().id_name() << "\n\n";
		++i;
	}
	return s;
}
