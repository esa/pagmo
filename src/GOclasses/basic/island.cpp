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

#include <boost/thread/thread.hpp>
#include <exception>

#include "../algorithms/go_algorithm.h"
#include "../problems/GOproblem.h"
#include "archipelago.h"
#include "island.h"
#include "population.h"

island::island(int n, GOProblem &p, const go_algorithm &al):m_pop(p,n),m_goa(al.clone()),m_a(0) {}

island::island(const island &i):m_pop(i.get_pop()),m_goa(i.m_goa->clone()),m_a(0) {}

island &island::operator=(const island &i)
{
	m_pop = i.get_pop();
	m_goa.reset(i.m_goa->clone());
	return *this;
}

island::~island()
{
	lock_type lock(m_mutex);
}

void island::evolve()
{
	boost::thread(evolver(this));
}

Population island::get_pop() const
{
	lock_type lock(m_mutex);
	return Population(m_pop);
}

void island::set_archipelago(archipelago *a)
{
	m_a = a;
}

void island::evolver::operator()()
{
	if (!m_i->m_mutex.try_lock()) {
		std::cout << "Cannot start evolution while already evolving.\n";
		return;
	}
	try {
		m_i->m_pop = m_i->m_goa->evolve(m_i->m_pop,m_i->m_pop.problem());
		std::cout << "Evolution finished, best fitness is: " << m_i->m_pop.extractBestIndividual().getFitness() << '\n';
	} catch (const std::exception &e) {
		std::cout << "Error during evolution: " << e.what() << '\n';
	}
	m_i->m_mutex.unlock();
}
