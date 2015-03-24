/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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

#ifndef PAGMO_BASE_ISLAND_H
#define PAGMO_BASE_ISLAND_H

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread/thread.hpp>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "config.h"
#include "algorithm/base.h"
#include "migration/base_r_policy.h"
#include "migration/base_s_policy.h"
#include "population.h"
#include "problem/base.h"
#include "serialization.h"
#include "types.h"

namespace pagmo
{
// Forward declarations.
class archipelago;
class base_island;

/// Alias for the shared pointer to a pagmo::base_island.
typedef boost::shared_ptr<base_island> base_island_ptr;

/// Base island class.
/**
 * This class incorporates a pagmo::population and a pagmo::algorithm::base used to evolve the population. Each time the evolve() (or evolve_t()) method is called,
 * a derived island class will execute the algorithm's evolve method on the population. The actual mechanism of launching the evolve method is defined in the
 * derived class - see \ref evolution_details "below" for more details. While evolution is undergoing, the island is locked down and no further operations will be allowed. The method join() can be used to wait until
 * evolution on the island has terminated. The busy() methods can be used to query the state of the island.
 *
 * The island can either exist as a stand-alone object or as a component of an archipelago.
 * If the island belongs to an archipelago, it can exchange individuals with other islands in the archipelago. The topology of the archipelago determines
 * the connections between islands, whereas every island can define migration policies to specify how to select and replace individuals during migration.
 * The relevant policy classes are migration::base_s_policy (selection policy) and migration::base_r_policy (replacement policy).
 * The probability of inserting migrating individuals into the island is regulated by the migration probability parameter. These migration parameters can be specified
 * upon island construction and they are given (hopefully) reasonable default values. See the constructors for detailed information.
 *
 * The interface of this class mirrors the interface of the population class. It is hence possible to get and set individuals, get the population size,
 * access the population champion, etc. The main difference
 * is that the methods of this class will never return references to internal members, in order to protect the internal state of the island while evolution is undergoing.
 * All getters methods will thus return copies instead of references, and all public methods will wait for an ongoing evolution to terminate before performing any action.
 *
 * \section evolution_details Implementation of the evolution methods
 *
 * When one of the evolution methods (evolve() or evolve_t()) is launched,
 * a local thread is opened and the perform_evolution() method is called from the new thread using as arguments the population and the algorithm stored in the island.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 * @author Marek Ruciński (marek.rucinski@gmail.com)
 */
class __PAGMO_VISIBLE base_island
{
	public:
		/// The archipelago class needs access to the internals of the island.
		friend class archipelago;
		/** @name Ctors, dtor and assignment operator.*/
		//@{
		base_island(const base_island &);
		explicit base_island(const algorithm::base &, const problem::base &, int,
			const migration::base_s_policy &,
			const migration::base_r_policy &);
		explicit base_island(const algorithm::base &, const population &,
			const migration::base_s_policy &,
			const migration::base_r_policy &);
		base_island &operator=(const base_island &);
		/// Clone method.
		/**
		 * Provided that the derived island implements properly the copy constructor, virtually all implementations of this method will
		 * look like this:
		 * \code
		 * return base_ptr(new derived_island(*this));
		 * \endcode
		 *
		 * @return pagmo::base_island_ptr to a copy of this.
		 */
		virtual base_island_ptr clone() const = 0;
		virtual ~base_island();
		//@}
		/** @name Input/output.*/
		//@{
		std::string human_readable_terse() const;
		std::string human_readable() const;
		virtual std::string get_name() const;
		//@}
		/** @name Evolution.
		 * Methods related to island evolution.
		 */
		//@{
		virtual void join() const;
		bool busy() const;
		void evolve(int = 1);
		void evolve_t(int);
		void interrupt();
		std::size_t get_evolution_time() const;
	protected:
		/// Method that implements the evolution of the population.
		virtual void perform_evolution(const algorithm::base &, population &) const = 0;
		virtual void thread_entry();
		virtual void thread_exit();
		//@}
	public:
		/** @name Getters and setters.*/
		//@{
		algorithm::base_ptr get_algorithm() const;
		void set_algorithm(const algorithm::base &);
		void set_x(population::size_type, const decision_vector &);
		void set_v(population::size_type, const decision_vector &);
		problem::base_ptr get_problem() const;
		population::size_type get_size() const;
		migration::base_s_policy_ptr get_s_policy() const;
		migration::base_r_policy_ptr get_r_policy() const;
		population get_population() const;
		void set_population(const population &);
		//@}
	private:
		// NOTE: in the next code line, it should be std::vector<std::pair<population::size_type, archipelago::size_type> >,
		// but this creates problems as at this point archipelago::siz_type is not defined and cannot be!!!
		std::vector<std::pair<population::size_type, population::size_type> > accept_immigrants(std::vector<std::pair<population::size_type, population::individual_type> > &);
		std::vector<population::individual_type> get_emigrants();
		// Evolver thread object. This is a callable helper object used to launch an evolution for a given number of iterations.
		struct int_evolver;
		// Time-dependent evolver thread object. This is a callable helper object used to launch an evolution for a specified amount of time.
		struct t_evolver;
		// RAII threads hook object.
		struct raii_thread_hook;
		friend struct raii_thread_hook;
	protected:
		/// Algorithm.
		algorithm::base_ptr			m_algo;
		/// Population.
		population				m_pop;
		/// Pointer that, if not null, points to the archipelago containing the island.
		archipelago				*m_archi;
		/// Total time spent by the island on evolution (in milliseconds).
		std::size_t				m_evo_time;
		/// Migration selection policy.
		migration::base_s_policy_ptr		m_s_policy;
		/// Migration replacement policy.
		migration::base_r_policy_ptr		m_r_policy;
		/// Evolution thread.
		boost::scoped_ptr<boost::thread>	m_evo_thread;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int version)
		{
			// Sync the island before doing anything.
			join();
			// TODO: Also, consider relation to save/load constructor data in island and mpi_island.
			ar & m_algo;
			ar & m_pop;
			ar & m_evo_time;
			ar & m_s_policy;
			ar & m_r_policy;
			boost::serialization::split_member(ar, *this, version);
		}
		template <class Archive>
		void save(Archive &, const unsigned int) const
		{}
		template <class Archive>
		void load(Archive &, const unsigned int)
		{
			// Upon loading we are going to set the archi pointer and the evo thread to 0.
			m_archi = 0;
			m_evo_thread.reset(0);
		}
};

std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const base_island &);

// Fake problem and algorithm used in the (de)serialization of island pointers.
namespace problem
{

class island_init: public base
{
	public:
		island_init():base(1) {}
		base_ptr clone() const
		{
			return base_ptr(new island_init(*this));
		}
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const {}
};

}

namespace algorithm
{

class island_init: public base
{
	public:
		island_init():base() {}
		base_ptr clone() const
		{
			return base_ptr(new island_init(*this));
		}
		void evolve(population &) const {}
};

}

}

BOOST_SERIALIZATION_ASSUME_ABSTRACT(pagmo::base_island)

#endif
