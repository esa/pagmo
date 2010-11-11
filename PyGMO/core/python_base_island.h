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

#ifndef PAGMO_PYTHON_BASE_ISLAND_H
#define PAGMO_PYTHON_BASE_ISLAND_H

#include <Python.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/python/class.hpp>
#include <boost/thread/thread.hpp>
#include <stdexcept>
#include <string>

#include "../../src/algorithm/base.h"
#include "../../src/base_island.h"
#include "../../src/config.h"
#include "../../src/exceptions.h"
#include "../../src/migration/base_r_policy.h"
#include "../../src/migration/base_s_policy.h"
#include "../../src/migration/best_s_policy.h"
#include "../../src/migration/fair_r_policy.h"
#include "../../src/population.h"
#include "../../src/problem/base.h"
#include "../../src/py_lock.h"
#include "../../src/serialization.h"
#include "../utils.h"

// Forward declarations.
namespace pagmo {

class python_base_island;

}

namespace boost { namespace serialization {

template <class Archive>
void save_construct_data(Archive &, const pagmo::python_base_island *, const unsigned int);

template <class Archive>
inline void load_construct_data(Archive &, pagmo::python_base_island *, const unsigned int);

}}

namespace pagmo {

// Base island class for re-implementation from Python.
class __PAGMO_VISIBLE python_base_island:  public base_island, public boost::python::wrapper<base_island>
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
		explicit python_base_island(const problem::base &prob, const algorithm::base &algo, int n = 0,
			const double &migr_prob = 1,
			const migration::base_s_policy &s_policy = migration::best_s_policy(),
			const migration::base_r_policy &r_policy = migration::fair_r_policy()):
			base_island(prob,algo,n,migr_prob,s_policy,r_policy) {}
		explicit python_base_island(const population &pop, const algorithm::base &algo,
			const double &migr_prob = 1,
			const migration::base_s_policy &s_policy = migration::best_s_policy(),
			const migration::base_r_policy &r_policy = migration::fair_r_policy()):
			base_island(pop,algo,migr_prob,s_policy,r_policy) {}
		// NOTE: why is this necessary?
		python_base_island(const base_island &isl):base_island(isl) {}
		~python_base_island()
		{
			// Call the re-implemented join().
			python_base_island::join();
		}
		python_base_island &operator=(const python_base_island &other)
		{
			base_island::operator=(other);
			return *this;
		}
		base_island_ptr clone() const
		{
			py_lock lock;
			base_island_ptr retval = this->get_override("__copy__")();
			if (!retval) {
				pagmo_throw(std::runtime_error,"island's __copy__() method returns a NULL pointer, please check the implementation");
			}
			return retval;
		}
		std::string get_name() const
		{
			py_lock lock;
			if (boost::python::override f = this->get_override("get_name")) {
				return f();
			}
			return base_island::get_name();
		}
		std::string default_get_name() const
		{
			return this->base_island::get_name();
		}
		// We need to reimplement this so that before attempting the thread
		// join we release the GIL. Otherwise, joining will block every other operation in the
		// separate threads calling from Python.
		void join() const
		{
			scoped_gil_release release;
			base_island::join();
		}
		population py_perform_evolution(algorithm::base_ptr a, const population &pop) const
		{
			py_lock lock;
			if (boost::python::override f = this->get_override("_perform_evolution")) {
				return f(a,pop);
			}
			pagmo_throw(not_implemented_error,"island's _perform_evolution method has not been implemented");
		}
	protected:
		// An island implemented in Python is never blocking: evolution goes into separate process.
		bool is_blocking_impl() const
		{
			return false;
		}
		void perform_evolution(const algorithm::base &a, population &pop) const
		{
			pop = py_perform_evolution(a.clone(),population(pop));
		}
	private:
		template <class Archive>
		friend void boost::serialization::save_construct_data(Archive &, const python_base_island *, const unsigned int);
		template <class Archive>
		friend void boost::serialization::load_construct_data(Archive &, python_base_island *, const unsigned int);
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_island>(*this);
			ar & boost::serialization::base_object<boost::python::wrapper<base_island> >(*this);
		}
};

struct python_base_island_pickle_suite : boost::python::pickle_suite
{
	static boost::python::tuple getinitargs(const python_base_island &isl)
	{
		pagmo::py_lock lock;
		return boost::python::make_tuple(isl.get_problem(),isl.get_algorithm());
	}
	static boost::python::tuple getstate(boost::python::object obj)
	{
		pagmo::py_lock lock;
		std::stringstream ss;
		const python_base_island &isl = boost::python::extract<python_base_island const &>(obj)();
		boost::archive::text_oarchive oa(ss);
		oa << isl;
		return boost::python::make_tuple(obj.attr("__dict__"),ss.str(),isl.get_algorithm());
	}
	static void setstate(boost::python::object obj, boost::python::tuple state)
	{
		pagmo::py_lock lock;
		if (len(state) != 3)
		{
			PyErr_SetObject(PyExc_ValueError,("expected 3-item tuple in call to __setstate__; got %s" % state).ptr());
			boost::python::throw_error_already_set();
		}
		// Restore the object's __dict__.
		boost::python::dict d = boost::python::extract<boost::python::dict>(obj.attr("__dict__"))();
		d.update(state[0]);
		// Restore the internal state of the C++ object.
		python_base_island &isl = boost::python::extract<python_base_island &>(obj)();
		const std::string str = boost::python::extract<std::string>(state[1]);
		std::stringstream ss(str);
		boost::archive::text_iarchive ia(ss);
		ia >> isl;
		const algorithm::base_ptr algo = boost::python::extract<algorithm::base_ptr>(state[2]);
		isl.set_algorithm(*algo);
	}
	static bool getstate_manages_dict()
	{
		return true;
	}
};

}

namespace boost { namespace serialization {

template <class Archive>
inline void save_construct_data(Archive &ar, const pagmo::python_base_island *isl, const unsigned int)
{
	// Save data required to construct instance.
	pagmo::problem::base_ptr prob = isl->m_pop.problem().clone();
	pagmo::algorithm::base_ptr algo = isl->m_algo->clone();
	ar << prob;
	ar << algo;
}

template <class Archive>
inline void load_construct_data(Archive &ar, pagmo::python_base_island *isl, const unsigned int)
{
	// Retrieve data from archive required to construct new instance.
	pagmo::problem::base_ptr prob;
	pagmo::algorithm::base_ptr algo;
	ar >> prob;
	ar >> algo;
	// Invoke inplace constructor to initialize instance of the island.
	::new(isl)pagmo::python_base_island(*prob,*algo);
}

}} //namespaces

BOOST_CLASS_EXPORT(pagmo::python_base_island);

#endif
