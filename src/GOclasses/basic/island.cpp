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

#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <boost/timer.hpp>
#include <stdexcept>

#include "../../exceptions.h"
#include "../algorithms/go_algorithm.h"
#include "../problems/GOproblem.h"
#include "archipelago.h"
#include "island.h"
#include "population.h"

size_t island::id_counter = 0;
boost::mutex island::id_mutex;

size_t island::get_new_id()
{
	lock_type lock(id_mutex);
	const size_t retval = id_counter;
	++id_counter;
	return retval;
}

island::island(int n, const GOProblem &p, const go_algorithm &al):m_id(get_new_id()),m_pop(p,n),m_goa(al.clone()),m_a(0) {}

island::island(const island &i):m_id(get_new_id()),m_pop(i.get_pop()),m_goa(i.m_goa->clone()),m_a(0) {}

island &island::operator=(const island &i)
{
	m_id = get_new_id();
	m_pop = i.get_pop();
	m_goa.reset(i.m_goa->clone());
	return *this;
}

island::~island()
{
	join();
}

void island::evolve(int N)
{
	if (!m_mutex.try_lock()) {
// WORKAROUND: apparently there are some issues here with  exception throwing
// under MinGW when evolution is running in another thread. This needs to be investigated.
#ifdef PAGMO_WIN32
		std::cout << "Cannot evolve while still evolving!\n";
		return;
#else
		pagmo_throw(std::runtime_error,"cannot evolve while still evolving");
#endif
	}
	boost::thread(int_evolver(this,N));
}

void island::evolve_t(const double &t)
{
	if (!m_mutex.try_lock()) {
#ifdef PAGMO_WIN32
		std::cout << "Cannot evolve while still evolving!\n";
		return;
#else
		pagmo_throw(std::runtime_error,"cannot evolve while still evolving");
#endif
	}
	boost::thread(t_evolver(this,t));
}

void island::join() const {
	lock_type lock(m_mutex);
}

bool island::active() const {
	if (!m_mutex.try_lock()) {
		return true;
	}
	m_mutex.unlock();
	return false;
}

size_t island::id() const
{
	return m_id;
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

void island::int_evolver::operator()()
{
	try {
		for (int i = 0; i < m_n; ++i) {
			m_i->m_pop = m_i->m_goa->evolve(m_i->m_pop);
			std::cout << "Evolution finished, best fitness is: " << m_i->m_pop.extractBestIndividual().getFitness() << '\n';
		}
	} catch (const std::exception &e) {
		std::cout << "Error during evolution: " << e.what() << '\n';
	} catch (...) {
		std::cout << "Unknown exception caught. :(\n";
	}
	m_i->m_mutex.unlock();
}

void island::t_evolver::operator()()
{
	try {
		boost::timer timer;
		while (true) {
			m_i->m_pop = m_i->m_goa->evolve(m_i->m_pop);
			std::cout << "Evolution finished, best fitness is: " << m_i->m_pop.extractBestIndividual().getFitness() << '\n';
			if (timer.elapsed() > m_t) {
				break;
			}
		}
	} catch (const std::exception &e) {
		std::cout << "Error during evolution: " << e.what() << '\n';
	} catch (...) {
		std::cout << "Unknown exception caught. :(\n";
	}
	m_i->m_mutex.unlock();
}
