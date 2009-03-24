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

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <iostream>
#include <stdexcept>

#include "../../atomic_counters/atomic_counters.h"
#include "../../exceptions.h"
#include "../algorithms/go_algorithm.h"
#include "../problems/GOproblem.h"
#include "archipelago.h"
#include "individual.h"
#include "island.h"
#include "population.h"

PaGMO::atomic_counter_size_t island::id_counter(1);

size_t island::get_new_id()
{
	return (size_t)(id_counter++);
}

island::island(const GOProblem &p, const go_algorithm &al)
		:m_id(get_new_id()),
		m_pop(p),
		m_goa(al.clone()),
		m_a(0),
		m_evo_time(0)
{
}

island::island(const GOProblem &p, const go_algorithm &al, int n)
		:m_id(get_new_id()),
		m_pop(p, n),
		m_goa(al.clone()),
		m_a(0),
		m_evo_time(0)
{
}

island::island(const GOProblem& p, const go_algorithm& al, int n, const MigrationSelectionPolicy& msp, const MigrationReplacementPolicy& mrp)
		:m_id(get_new_id()),
		m_pop(p, n),
		m_goa(al.clone()),
		m_a(0),
		m_evo_time(0),
		migrationSelectionPolicy(msp.clone()),
		migrationReplacementPolicy(mrp.clone())
{
}

island::island(const island &i)
		:m_id(get_new_id()),
		m_pop(i.population()),
		m_goa(i.m_goa->clone()),
		m_a(0),
		m_evo_time(i.m_evo_time),
		migrationSelectionPolicy(i.migrationSelectionPolicy ? i.migrationSelectionPolicy->clone() : 0),
		migrationReplacementPolicy(i.migrationReplacementPolicy ? i.migrationReplacementPolicy->clone() : 0)
{
}

island &island::operator=(const island &i)
{
	if (this != &i) {
		m_pop = i.population();
		m_goa.reset(i.m_goa->clone());
		m_evo_time = i.m_evo_time;
		//TODO: policies!!!
	}
	return *this;
}

island::~island()
{
	join();
}


Population island::population() const
{
	join();
	return m_pop;
}

const GOProblem &island::problem() const
{
	join();
	return m_pop.problem();
}

const go_algorithm &island::algorithm() const
{
	join();
	return *m_goa;
}

void island::set_algorithm(const go_algorithm &a)
{
	join();
	m_goa.reset(a.clone());
}

const MigrationSelectionPolicy* island::getMigrationSelectionPolicy() const
{
	join();
	return migrationSelectionPolicy.get();
}

void island::setMigrationSelectionPolicy(const MigrationSelectionPolicy* msp)
{
	join();
	migrationSelectionPolicy.reset(msp ? msp->clone() : 0);
}

const MigrationReplacementPolicy* island::getMigrationReplacementPolicy() const
{
	join();
	return migrationReplacementPolicy.get();
}

void island::setMigrationReplacementPolicy(const MigrationReplacementPolicy* mrp)
{
	join();
	migrationReplacementPolicy.reset(mrp ? mrp->clone() : 0);
}

size_t island::id() const
{
	return m_id;
}

size_t island::size() const
{
	join();
	return m_pop.size();
}

size_t island::evo_time() const
{
	join();
	return m_evo_time;
}


Individual island::operator[](int n) const
{
	join();
	return m_pop[n];
}

void island::set_individual(int n, const Individual &i)
{
	join();
	m_pop.setIndividual(n,i);
}

void island::push_back(const Individual &i)
{
	join();
	m_pop.push_back(i);
}

void island::insert(int n, const Individual &i)
{
	join();
	m_pop.insert(n,i);
}

void island::erase(int n)
{
	join();
	m_pop.erase(n);
}

double island::mean() const
{
	join();
	return m_pop.evaluateMean();
}

double island::std() const
{
	join();
	return m_pop.evaluateStd();
}

Individual island::best() const
{
	join();
	return m_pop.extractBestIndividual();
}

Individual island::worst() const
{
	join();
	return m_pop.extractBestIndividual();
}



std::vector<Individual> island::getMigratingIndividuals()
{
	if(migrationSelectionPolicy) {
		return migrationSelectionPolicy->selectForMigration(m_pop);
	} else {
		/// \todo Throw an exception here?
		//return an empty population
		return std::vector<Individual>();
	}
}


void island::acceptMigratingIndividuals(const std::vector<Individual>& incomingPopulation)
{
	if(migrationReplacementPolicy) {
		std::list<std::pair<int, int> > replacement = migrationReplacementPolicy->selectForReplacement(incomingPopulation, m_pop);
		
		for(std::list<std::pair<int, int> >::iterator rep = replacement.begin(); rep != replacement.end(); rep++) {
			/// \todo This method does some strange checking - maybe we could ommit it?
			m_pop.setIndividual((*rep).first, incomingPopulation[(*rep).second]);
		}
	} else {
		/// \todo Throw an exception here?
	}	
}



void island::evolve(int N)
{
	if (m_evo_mutex.try_lock()) {
		try {
			boost::thread(int_evolver(this,N));
		} catch(...) {
			std::cout << "Failed to launch the thread: " << m_id << std::endl;
			/// \todo throw;
		}
	} else {
// WORKAROUND: apparently there are some issues here with  exception throwing
// under MinGW when evolution is running in another thread. This needs to be investigated.
#ifdef PAGMO_WIN32
		std::cout << "Cannot evolve while still evolving!\n";
#else
		pagmo_throw(std::runtime_error,"cannot evolve while still evolving");
#endif
	}
}

void island::evolve_t(const size_t &t)
{
	if (m_evo_mutex.try_lock()) {
		try {
			boost::thread(t_evolver(this,t));
		} catch(...) {
			std::cout << "Failed to launch the thread: " << m_id << std::endl;
			/// \todo throw;
		}
	} else {
#ifdef PAGMO_WIN32
		std::cout << "Cannot evolve while still evolving!\n";
#else
		pagmo_throw(std::runtime_error,"cannot evolve while still evolving");
#endif
	}
}

void island::join() const {
	lock_type lock(m_evo_mutex);
}

bool island::busy() const {
	if (!m_evo_mutex.try_lock()) {
		return true;
	}
	m_evo_mutex.unlock();
	return false;
}


void island::set_archipelago(archipelago *a)
{
	m_a = a;
}

void island::int_evolver::operator()()
{
	const boost::posix_time::ptime start = boost::posix_time::microsec_clock::local_time();
	try {
		for (int i = 0; i < m_n; ++i) {
			if (m_i->m_a) {
				//lock_type lock(m_i->m_topo_mutex);
				m_i->m_a->preEvolutionCallback(*m_i);
				m_i->m_pop.problem().pre_evolution();
			}
			m_i->m_pop = m_i->m_goa->evolve(m_i->m_pop);
			//std::cout << "Evolution finished, best fitness is: " << m_i->m_pop.extractBestIndividual().getFitness() << '\n';
			if (m_i->m_a) {
				//lock_type lock(m_i->m_topo_mutex);
				m_i->m_a->postEvolutionCallback(*m_i);
				m_i->m_pop.problem().post_evolution();
			}
		}
	} catch (const std::exception &e) {
		std::cout << "Error during evolution: " << e.what() << '\n';
	} catch (...) {
		std::cout << "Unknown exception caught. :(\n";
	}
	// We must take care of potentially low-accuracy clocks, where the time difference could be negative for
	// _really_ short evolution times. In that case do not add anything to the total evolution time.
	const boost::posix_time::time_duration diff = boost::posix_time::microsec_clock::local_time() - start;
	if (diff.total_milliseconds() >= 0) {
		m_i->m_evo_time += diff.total_milliseconds();
	}
	m_i->m_evo_mutex.unlock();
}

// Perform at least one evolution, and continue evolving until at least a certain amount of time has passed.
void island::t_evolver::operator()()
{
	const boost::posix_time::ptime start = boost::posix_time::microsec_clock::local_time();
	boost::posix_time::time_duration diff;
	try {
		do {
			if (m_i->m_a) {
				//lock_type lock(m_i->m_topo_mutex);
				//m_i->m_a->m_top->pre_evolution(*m_i); //Call migration scheme
				m_i->m_pop.problem().pre_evolution();
			}
			m_i->m_pop = m_i->m_goa->evolve(m_i->m_pop);
			diff = boost::posix_time::microsec_clock::local_time() - start;
			//std::cout << "Evolution finished, best fitness is: " << m_i->m_pop.extractBestIndividual().getFitness() << '\n';
			if (m_i->m_a) {
				//lock_type lock(m_i->m_topo_mutex);
				//m_i->m_a->m_top->post_evolution(*m_i); //Call migration scheme
				m_i->m_pop.problem().post_evolution();
			}
			// Take care of negative timings.
		} while (diff.total_milliseconds() < 0 || (size_t)diff.total_milliseconds() < m_t);
	} catch (const std::exception &e) {
		std::cout << "Error during evolution: " << e.what() << '\n';
	} catch (...) {
		std::cout << "Unknown exception caught. :(\n";
	}
	m_i->m_evo_time += diff.total_milliseconds();
	m_i->m_evo_mutex.unlock();
}


std::ostream &operator<<(std::ostream &s, const island &isl)
{
	s << "ID:                 " << isl.id() << '\n';
	s << "Population size:    " << isl.size() << '\n';
	s << "Evolution time:     " << isl.evo_time() << '\n';
	s << "Algorithm type:     " << isl.algorithm().id_name() << '\n';
	s << "Selection policy:   ";
	if(isl.migrationSelectionPolicy) {
		s << std::endl << *(isl.migrationSelectionPolicy);
	} else {
		s << "none" << std::endl;		
	}
	s << "Replacement policy: ";
	if(isl.migrationReplacementPolicy) {
		s << std::endl << *(isl.migrationReplacementPolicy);
	} else {
		s << "none" << std::endl;		
	}	
	boost::lock_guard<boost::mutex> lock(isl.m_evo_mutex);
	s << isl.m_pop;
	return s;
}
