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

#ifndef PAGMO_PYTHON_BASE_ISLAND_H
#define PAGMO_PYTHON_BASE_ISLAND_H

#include <Python.h>
#include <boost/python/class.hpp> // For pickle suite.
#include <boost/python/dict.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/object.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/ptr.hpp>
#include <boost/python/wrapper.hpp>
#include <stdexcept>
#include <sstream>
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
#include "../../src/serialization.h"
#include "../../src/rng.h"
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
		// RAII gil releaser. See:
		// http://wiki.python.org/moin/boost.python/HowTo#MultithreadingSupportformyfunction
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
		explicit python_base_island(const algorithm::base &algo, const problem::base &prob, int n = 0,
			const migration::base_s_policy &s_policy = migration::best_s_policy(),
			const migration::base_r_policy &r_policy = migration::fair_r_policy()):
			base_island(algo,prob,n,s_policy,r_policy),
			boost::python::wrapper<base_island>(),
			m_gstate() {}
		explicit python_base_island(const algorithm::base &algo, const population &pop,
			const migration::base_s_policy &s_policy = migration::best_s_policy(),
			const migration::base_r_policy &r_policy = migration::fair_r_policy()):
			base_island(algo,pop,s_policy,r_policy),
			boost::python::wrapper<base_island>(),
			m_gstate() {}
		python_base_island(const python_base_island &isl):base_island(isl), boost::python::wrapper<base_island>(), m_gstate() {}
		~python_base_island()
		{
			// Call the re-implemented join().
			python_base_island::join();
			pagmo_assert(m_gstate == PyGILState_STATE());
		}
		python_base_island &operator=(const python_base_island &other)
		{
			base_island::operator=(other);
			pagmo_assert(m_gstate == PyGILState_STATE());
			return *this;
		}
		base_island_ptr clone() const
		{
			join();
			base_island_ptr retval = this->get_override("__get_deepcopy__")();
			if (!retval) {
				pagmo_throw(std::runtime_error,"island's __get_deepcopy__() method returns a NULL pointer, please check the implementation");
			}
			return retval;
		}
		std::string get_name() const
		{
			join();
			if (boost::python::override f = this->get_override("get_name")) {
			#if BOOST_WORKAROUND(BOOST_MSVC, <= 1700)
				return boost::python::call<std::string>(this->get_override("get_name").ptr());
			#else
				return f();				
			#endif
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
		population py_perform_evolution(const algorithm::base *a, const population &pop) const
		{
			if (boost::python::override f = this->get_override("_perform_evolution")) {
				return f(boost::python::ptr(a),pop);
			}
			pagmo_throw(not_implemented_error,"island's _perform_evolution method has not been implemented");
		}
	protected:
		void perform_evolution(const algorithm::base &a, population &pop) const
		{
			// note that we pass 'the real' a 
			population retval(py_perform_evolution(&a,pop));
			a.reset_rngs(rng_generator::get<rng_uint32>()());
			// Check that the implementation of the evolve method in Python did not screw up the problem.
			if (pop.problem() != retval.problem()) {
				pagmo_throw(std::runtime_error,"the island's perform_evolution method returned a population whose problem is inconsistent with that of the input population");
			}
			pop = retval;
		}
		void thread_entry()
		{
			// PyGILState_* functions from PYthon >= 2.3. See:
			// http://docs.python.org/c-api/init.html
			m_gstate = PyGILState_Ensure();
		}
		void thread_exit()
		{
			PyGILState_Release(m_gstate);
			m_gstate = PyGILState_STATE();
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
		// The only data member is the GIL state variable.
		PyGILState_STATE m_gstate;
};

struct python_base_island_pickle_suite : boost::python::pickle_suite
{
	static boost::python::tuple getinitargs(const python_base_island &isl)
	{
		return boost::python::make_tuple(isl.get_algorithm(),isl.get_problem());
	}
	static boost::python::tuple getstate(boost::python::object obj)
	{
		const python_base_island &isl = boost::python::extract<python_base_island const &>(obj)();
		std::stringstream ss;
		boost::archive::text_oarchive oa(ss);
		oa << isl;
		return boost::python::make_tuple(obj.attr("__dict__"),ss.str(),isl.get_algorithm(),isl.get_population());
	}
	static void setstate(boost::python::object obj, boost::python::tuple state)
	{
		if (len(state) != 4)
		{
			PyErr_SetObject(PyExc_ValueError,("expected 4-item tuple in call to __setstate__; got %s" % state).ptr());
			boost::python::throw_error_already_set();
		}
		python_base_island &isl = boost::python::extract<python_base_island &>(obj)();
		// Restore the object's __dict__.
		boost::python::dict d = boost::python::extract<boost::python::dict>(obj.attr("__dict__"))();
		d.update(state[0]);
		// Restore the internal state of the C++ object.
		const std::string str = boost::python::extract<std::string>(state[1]);
		std::stringstream ss(str);
		boost::archive::text_iarchive ia(ss);
		ia >> isl;
		// Restore separately the algorithm and the population.
		// NOTE: here (and elsewhere in similar situations) we could avoid the need to deal separately with population and/or algorithm:
		// as long as we are not dealing with Python-extended objects, we are sure that C++ serialization is enough. Optimize like this
		// in case serialization eventually turns out to be a bottleneck.
		const algorithm::base_ptr algo = boost::python::extract<algorithm::base_ptr>(state[2]);
		isl.set_algorithm(*algo);
		const population pop = boost::python::extract<population>(state[3]);
		isl.set_population(pop);
	}
	static bool getstate_manages_dict()
	{
		return true;
	}
};

}

namespace boost { namespace serialization {

// Do no need to save any data, will use fake problem and algorithm for pointer initialization.
template <class Archive>
inline void save_construct_data(Archive &, const pagmo::python_base_island *, const unsigned int)
{}

template <class Archive>
inline void load_construct_data(Archive &, pagmo::python_base_island *isl, const unsigned int)
{
	// Invoke inplace constructor to initialize instance of the island.
	::new(isl)pagmo::python_base_island(pagmo::algorithm::island_init(),pagmo::problem::island_init());
}

}} //namespaces

BOOST_CLASS_EXPORT(pagmo::python_base_island)

#endif
