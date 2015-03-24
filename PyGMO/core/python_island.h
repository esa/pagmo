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
 *   the Free Software Foundation; either version 3 of the License, or       *
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

#ifndef PAGMO_PYTHON_ISLAND_H
#define PAGMO_PYTHON_ISLAND_H

#include <Python.h>
#include <typeinfo>
#include "../../src/config.h"
#include "../../src/base_island.h"
#include "../../src/island.h"
#include "../../src/algorithm/base.h"
#include "../../src/migration/base_r_policy.h"
#include "../../src/migration/base_s_policy.h"
#include "../../src/migration/best_s_policy.h"
#include "../../src/migration/fair_r_policy.h"
#include "../../src/population.h"
#include "../../src/problem/base.h"
#include "../../src/serialization.h"
#include "../algorithm/python_base.h"
#include "../problem/python_base.h"

// Forward declarations.
namespace pagmo {

class python_island;

}

namespace boost { namespace serialization {

template <class Archive>
void save_construct_data(Archive &, const pagmo::python_island *, const unsigned int);

template <class Archive>
inline void load_construct_data(Archive &, pagmo::python_island *, const unsigned int);

}}

namespace pagmo {

// We need to reimplement the local island class to handle the case in which the island is instantiated with
// algorithms/problems implemented in Python and to re-implement the join() method not to block other Python
// computations.
class __PAGMO_VISIBLE python_island: public island
{
		// RAII gil releaser.
		class scoped_gil_release
		{
			public:
				scoped_gil_release()
				{
					m_thread_state = PyEval_SaveThread();
				}
				~scoped_gil_release()
				{
					PyEval_RestoreThread(m_thread_state);
					m_thread_state = NULL;
				}
			private:
				PyThreadState *m_thread_state;
		};
	public:
		explicit python_island(const algorithm::base &algo, const problem::base &prob, int n = 0,
			const migration::base_s_policy &s_policy = migration::best_s_policy(),
			const migration::base_r_policy &r_policy = migration::fair_r_policy()):
			island(algo,prob,n,s_policy,r_policy),m_gstate() {}
		explicit python_island(const algorithm::base &algo, const population &pop,
			const migration::base_s_policy &s_policy = migration::best_s_policy(),
			const migration::base_r_policy &r_policy = migration::fair_r_policy()):
			island(algo,pop,s_policy,r_policy),m_gstate() {}
		python_island(const python_island &isl):island(isl),m_gstate() {}
		~python_island()
		{
			// Call the re-implemented join().
			python_island::join();
			pagmo_assert(m_gstate == PyGILState_STATE());
		}
		python_island &operator=(const python_island &other)
		{
			island::operator=(other);
			pagmo_assert(m_gstate == PyGILState_STATE());
			return *this;
		}
		base_island_ptr clone() const
		{
			return base_island_ptr(new python_island(*this));
		}
		void join() const
		{
			scoped_gil_release release;
			base_island::join();
		}
	protected:
		void thread_entry()
		{
			if (is_pythonic()) {
				m_gstate = PyGILState_Ensure();
			}
		}
		void thread_exit()
		{
			if (is_pythonic()) {
				PyGILState_Release(m_gstate);
				m_gstate = PyGILState_STATE();
			}
		}
	public:
		bool is_pythonic() const
		{
			int n_pythonic_items = 0;
			try {
				algorithm::python_base trial = dynamic_cast<algorithm::python_base &>(*m_algo);
				++n_pythonic_items;
			} catch (const std::bad_cast &) {}
			try {
				const problem::python_base trial = dynamic_cast<problem::python_base const &>(m_pop.problem());
				++n_pythonic_items;
			} catch (const std::bad_cast &) {}
			return (n_pythonic_items > 0);
		}
	private:
		template <class Archive>
		friend void boost::serialization::save_construct_data(Archive &, const python_island *, const unsigned int);
		template <class Archive>
		friend void boost::serialization::load_construct_data(Archive &, python_island *, const unsigned int);
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			// Join will be done here already.
			ar & boost::serialization::base_object<island>(*this);
		}
		// The only data member is the GIL state variable.
		PyGILState_STATE m_gstate;
};

}

namespace boost { namespace serialization {

template <class Archive>
inline void save_construct_data(Archive &, const pagmo::python_island *, const unsigned int)
{}

template <class Archive>
inline void load_construct_data(Archive &, pagmo::python_island *isl, const unsigned int)
{
	::new(isl)pagmo::python_island(pagmo::algorithm::island_init(),pagmo::problem::island_init());
}

}} //namespaces

BOOST_CLASS_EXPORT(pagmo::python_island)

#endif
