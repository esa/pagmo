#ifndef PAGMO_ATOMIC_COUNTER_GENERIC_H
#define PAGMO_ATOMIC_COUNTER_GENERIC_H

#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>

#include "base_atomic_counter.h"

namespace PaGMO
{
	template <class IntType>
	class atomic_counter_generic: public base_atomic_counter<IntType,atomic_counter_generic<IntType> >
	{
			typedef base_atomic_counter<IntType,atomic_counter_generic<IntType> > ancestor;
		public:
			atomic_counter_generic():ancestor::base_atomic_counter(),m_mutex() {}
			template <class IntType2>
			atomic_counter_generic &operator+=(const IntType2 &n) {
				boost::lock_guard<boost::mutex> lock(m_mutex);
				this->m_value += n;
				return *this;
			}
			template <class IntType2>
			atomic_counter_generic &operator-=(const IntType2 &n) {
				boost::lock_guard<boost::mutex> lock(m_mutex);
				this->m_value -= n;
				return *this;
			}
			template <class IntType2>
			bool compare_and_swap(const IntType2 &oldval, const IntType2 &newval) {
				boost::lock_guard<boost::mutex> lock(m_mutex);
				if (this->m_value == oldval) {
					this->m_value = newval;
					return true;
				}
				return false;
			}
		private:
			boost::mutex m_mutex;
	};
}

#endif
