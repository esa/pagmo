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

// 04/01/2009: Initial version by Francesco Biscani.

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <iostream>
#include <stdexcept>

#include "../../atomic_counters/atomic_counters.h"
#include "../../exceptions.h"
#include "../algorithms/base.h"
#include "../problems/base.h"
#include "archipelago.h"
#include "individual.h"
#include "island.h"
#include "population.h"

namespace pagmo
{

atomic_counter_size_t island::id_counter(1);

size_t island::get_new_id()
{
	return (id_counter++).get_value();
}

island::island(const problem::base &p, const algorithm::base &al)
		:m_id(get_new_id()),
		m_pop(p),
		m_goa(al.clone()),
		m_a(0),
		m_evo_time(0)
{
}

island::island(const problem::base &p, const algorithm::base &al, int n)
		:m_id(get_new_id()),
		m_pop(p, n),
		m_goa(al.clone()),
		m_a(0),
		m_evo_time(0)
{
}

island::island(const problem::base& p, const algorithm::base& al, int n, const MigrationPolicy& mp)
		:m_id(get_new_id()),
		m_pop(p, n),
		m_goa(al.clone()),
		m_a(0),
		m_evo_time(0),
		migrationPolicy(new MigrationPolicy(mp))
{
}

island::island(const island &i)
		:m_id(get_new_id()),
		m_pop(i.get_population()),
		m_goa(i.m_goa->clone()),
		m_a(0),
		m_evo_time(i.m_evo_time),
		migrationPolicy(i.migrationPolicy ? new MigrationPolicy(*(i.migrationPolicy)) : 0)
{
}

island &island::operator=(const island &i)
{
	if (this != &i) {
		m_pop = i.get_population();
		m_goa.reset(i.m_goa->clone());
		m_evo_time = i.m_evo_time;
		/** TODO Make it consistent with the copy constructor!!! */
	}
	return *this;
}

island::~island()
{
	join();
}


population island::get_population() const
{
	join();
	return m_pop;
}

const problem::base &island::problem() const
{
	join();
	return m_pop.problem();
}

const algorithm::base &island::algorithm() const
{
	join();
	return *m_goa;
}

void island::set_algorithm(const algorithm::base &a)
{
	join();
	m_goa.reset(a.clone());
}

const MigrationSelectionPolicy& island::getMigrationSelectionPolicy() const
{
	join();
	if (!migrationPolicy) {
		pagmo_throw(value_error, "The island has no associated migration policy!");
	}
	return migrationPolicy->getMigrationSelectionPolicy();
}

void island::setMigrationSelectionPolicy(const MigrationSelectionPolicy& msp)
{
	join();
	if (!migrationPolicy) {
		pagmo_throw(value_error, "The island has no associated migration policy!");
	}
	migrationPolicy->setMigrationSelectionPolicy(msp);
}

const MigrationReplacementPolicy& island::getMigrationReplacementPolicy() const
{
	join();
	if (!migrationPolicy) {
		pagmo_throw(value_error, "The island has no associated migration policy!");
	}
	return migrationPolicy->getMigrationReplacementPolicy();
}

void island::setMigrationReplacementPolicy(const MigrationReplacementPolicy& mrp)
{
	join();
	if (!migrationPolicy) {
		pagmo_throw(value_error, "The island has no associated migration policy!");
	}
	migrationPolicy->setMigrationReplacementPolicy(mrp);
}

const MigrationPolicy& island::getMigrationPolicy() const
{
	join();
	if (!migrationPolicy) {
		pagmo_throw(value_error, "The island has no associated migration policy!");
	}
	return *migrationPolicy;
}

void island::setMigrationPolicy(const MigrationPolicy* mp)
{
	join();
	migrationPolicy.reset(mp ? new MigrationPolicy(*mp) : 0);
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


individual island::operator[](int n) const
{
	join();
	return m_pop[n];
}

void island::set_individual(int n, const individual &i)
{
	join();
	m_pop.setIndividual(n,i);
}

void island::push_back(const individual &i)
{
	join();
	m_pop.push_back(i);
}

void island::insert(int n, const individual &i)
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

individual island::best() const
{
	join();
	return m_pop.extractBestIndividual();
}

individual island::worst() const
{
	join();
	return m_pop.extractWorstIndividual();
}


std::vector<individual> island::getMigratingIndividuals()
{
	if (migrationPolicy) {
		return migrationPolicy->getMigrationSelectionPolicy().selectForMigration(m_pop);
	} else {
		/// \todo Throw an exception here?
		//return an empty population
		return std::vector<individual>();
	}
}


void island::acceptMigratingIndividuals(const std::vector<individual>& incomingPopulation)
{
	if (migrationPolicy) {
		std::list<std::pair<int, int> > replacement = migrationPolicy->getMigrationReplacementPolicy().selectForReplacement(incomingPopulation, m_pop);

		for (std::list<std::pair<int, int> >::iterator rep = replacement.begin(); rep != replacement.end(); rep++) {
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
		} catch (...) {
			std::cout << "Failed to launch the thread: " << m_id << std::endl;
			/// \todo throw;
		}
	} else {
		pagmo_throw(std::runtime_error,"cannot evolve while still evolving");
	}
}

void island::evolve_t(const size_t &t)
{
	if (m_evo_mutex.try_lock()) {
		try {
			boost::thread(t_evolver(this,t));
		} catch (...) {
			std::cout << "Failed to launch the thread: " << m_id << std::endl;
			/// \todo throw;
		}
	} else {
		pagmo_throw(std::runtime_error,"cannot evolve while still evolving");
	}
}

void island::join() const
{
	lock_type lock(m_evo_mutex);
}

bool island::busy() const
{
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
		//Synchronise start with all other threads
		if (m_i->m_a) {
			m_i->m_a->sync_island_start();
		}

		for (int i = 0; i < m_n; ++i) {
			if (m_i->m_a) {
				//lock_type lock(m_i->m_topo_mutex);
				m_i->m_a->preEvolutionCallback(*m_i);
				m_i->m_pop.problem().pre_evolution(m_i->m_pop);
			}
			m_i->m_pop = m_i->m_goa->evolve(m_i->m_pop);
			//std::cout << "Evolution finished, best fitness is: " << m_i->m_pop.extractBestindividual().get_fitness() << '\n';
			if (m_i->m_a) {
				//lock_type lock(m_i->m_topo_mutex);
				m_i->m_a->postEvolutionCallback(*m_i);
				m_i->m_pop.problem().post_evolution(m_i->m_pop);
			}
		}
	} catch (const std::exception &e) {
		std::cout << "Error during island evolution: " << e.what() << '\n';
	} catch (...) {
		std::cout << "Error during island evolution, unknown exception caught. :(\n";
	}
	// We must take care of potentially low-accuracy clocks, where the time difference could be negative for
	// _really_ short evolution times. In that case do not add anything to the total evolution time.
	const boost::posix_time::time_duration diff = boost::posix_time::microsec_clock::local_time() - start;
	if (diff.total_milliseconds() >= 0) {
		m_i->m_evo_time += (size_t)(diff.total_milliseconds());
	}
	m_i->m_evo_mutex.unlock();
}

// Perform at least one evolution, and continue evolving until at least a certain amount of time has passed.
void island::t_evolver::operator()()
{
	const boost::posix_time::ptime start = boost::posix_time::microsec_clock::local_time();
	boost::posix_time::time_duration diff;
	try {
		// Synchronise start
		if (m_i->m_a) {
			m_i->m_a->sync_island_start();
		}

		do {
			if (m_i->m_a) {
				//lock_type lock(m_i->m_topo_mutex);
				//m_i->m_a->m_top->pre_evolution(*m_i); //Call migration scheme
				m_i->m_pop.problem().pre_evolution(m_i->m_pop);
			}
			m_i->m_pop = m_i->m_goa->evolve(m_i->m_pop);
			diff = boost::posix_time::microsec_clock::local_time() - start;
			//std::cout << "Evolution finished, best fitness is: " << m_i->m_pop.extractBestindividual().get_fitness() << '\n';
			if (m_i->m_a) {
				//lock_type lock(m_i->m_topo_mutex);
				//m_i->m_a->m_top->post_evolution(*m_i); //Call migration scheme
				m_i->m_pop.problem().post_evolution(m_i->m_pop);
			}
			// Take care of negative timings.
		} while (diff.total_milliseconds() < 0 || (size_t)diff.total_milliseconds() < m_t);
	} catch (const std::exception &e) {
		std::cout << "Error during evolution: " << e.what() << '\n';
	} catch (...) {
		std::cout << "Unknown exception caught. :(\n";
	}
	m_i->m_evo_time += (size_t)(diff.total_milliseconds());
	m_i->m_evo_mutex.unlock();
}


std::ostream &operator<<(std::ostream &s, const island &isl)
{
	s << "ID:                          " << isl.id() << '\n';
	s << "Population size:             " << isl.size() << '\n';
	s << "Evolution time:              " << isl.evo_time() << '\n';
	s << "Algorithm type:              " << isl.algorithm().id_name() << '\n';
	s << "Migration policy:            ";
	if (isl.migrationPolicy) {
		s << std::endl << *(isl.migrationPolicy);
	} else {
		s << "none" << std::endl;
	}
	boost::lock_guard<boost::mutex> lock(isl.m_evo_mutex);
	s << isl.m_pop;
	return s;
}

}
