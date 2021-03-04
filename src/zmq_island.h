/*****************************************************************************
 *	 Copyright (C) 2004-2015 The PaGMO development team,					 *
 *	 Advanced Concepts Team (ACT), European Space Agency (ESA)				 *
 *																			 *
 *	 https://github.com/esa/pagmo											 *
 *																			 *
 *	 act@esa.int															 *
 *																			 *
 *	 This program is free software; you can redistribute it and/or modify	 *
 *	 it under the terms of the GNU General Public License as published by	 *
 *	 the Free Software Foundation; either version 2 of the License, or		 *
 *	 (at your option) any later version.									 *
 *																			 *
 *	 This program is distributed in the hope that it will be useful,		 *
 *	 but WITHOUT ANY WARRANTY; without even the implied warranty of			 *
 *	 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the			 *
 *	 GNU General Public License for more details.							 *
 *																			 *
 *	 You should have received a copy of the GNU General Public License		 *
 *	 along with this program; if not, write to the							 *
 *	 Free Software Foundation, Inc.,										 *
 *	 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.				 *
 *****************************************************************************/

#ifndef PAGMO_ZMQ_ISLAND_H
#define PAGMO_ZMQ_ISLAND_H

#include <boost/scoped_ptr.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/thread/mutex.hpp>
#include <redox.hpp>
#include <zmq.hpp>
#include <list>
#include <set>
#include <string>

#include "base_island.h"
#include "config.h"
#include "algorithm/base.h"
#include "migration/base_r_policy.h"
#include "migration/base_s_policy.h"
#include "migration/best_s_policy.h"
#include "migration/fair_r_policy.h"
#include "population.h"
#include "problem/base.h"
#include "serialization.h"

// Forward declarations.
namespace pagmo {

class zmq_island;

}

namespace boost { namespace serialization {

template <class Archive>
void save_construct_data(Archive &, const pagmo::zmq_island *, const unsigned int);

template <class Archive>
inline void load_construct_data(Archive &, pagmo::zmq_island *, const unsigned int);

}}

namespace pagmo
{

/// ZMQ island class.
/**
 * This island can communicate with other ZeroMQ islands active on as well as local islands.
 * The intended use of this island class is to exchange solutions with remote archipelagos
 * through a ZeroMQ socket. The functionality is similar to the MPI island (TODO: add ref),
 * but using a different transport protocol.
 *
 * <b>NOTE</b>: this class is available only if PaGMO was compiled with ZeroMQ support.
 *
 * @author Jose Diez (me@jdiez.me)
 */
class __PAGMO_VISIBLE zmq_island: public base_island
{
		template <class Archive>
		friend void boost::serialization::save_construct_data(Archive &, const pagmo::zmq_island *, const unsigned int);
		template <class Archive>
		friend void boost::serialization::load_construct_data(Archive &, pagmo::zmq_island *, const unsigned int);
	public:
		zmq_island(const zmq_island &);
	    ~zmq_island();
		explicit zmq_island(const algorithm::base &, const problem::base &, int = 0,
			const migration::base_s_policy & = migration::best_s_policy(),
			const migration::base_r_policy & = migration::fair_r_policy());
		explicit zmq_island(const algorithm::base &, const population &,
			const migration::base_s_policy & = migration::best_s_policy(),
			const migration::base_r_policy & = migration::fair_r_policy());
		zmq_island &operator=(const zmq_island &);
		base_island_ptr clone() const;
	protected:
		void perform_evolution(const algorithm::base &, population &) const;
	public:
		typedef void (*callback)(zmq::message_t&);

		std::string get_name() const;

		void set_broker_details(std::string, int);
		void set_token(std::string);
		bool connect();
		void set_evolve(bool);
		bool get_evolve();

		void set_callback(zmq_island::callback);
		void disable_callback();
		void set_ip(std::string ip);
		void close();

	protected:
		std::string m_brokerHost;
		int			m_brokerPort;
		std::string m_token;
		std::string m_IP;
		int			m_localPort;

		bool m_initialised;
		bool m_evolve;

		zmq_island::callback m_callback;

	private:
		redox::Redox m_brokerConn;
		redox::Subscriber m_brokerSubscriber;

		zmq::context_t m_zmqContext;
		zmq::socket_t  m_publisherSocket;
		zmq::socket_t  m_subscriptionSocket;

		void connect_host(std::string);
		void disconnect();

		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			// Join is already done in base_island.
			ar & boost::serialization::base_object<base_island>(*this);
		}
};

}

namespace boost { namespace serialization {

template <class Archive>
inline void save_construct_data(Archive &ar, const pagmo::zmq_island *isl, const unsigned int)
{
	// Save data required to construct instance.
	pagmo::algorithm::base_ptr algo = isl->m_algo->clone();
	pagmo::problem::base_ptr prob = isl->m_pop.problem().clone();
	ar << algo;
	ar << prob;
}

template <class Archive>
inline void load_construct_data(Archive &ar, pagmo::zmq_island *isl, const unsigned int)
{
	// Retrieve data from archive required to construct new instance.
	pagmo::algorithm::base_ptr algo;
	pagmo::problem::base_ptr prob;
	ar >> algo;
	ar >> prob;
	// Invoke inplace constructor to initialize instance of the algorithm.
	::new(isl)pagmo::zmq_island(*algo,*prob);
}

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::zmq_island)

#endif
