#ifndef PAGMO_ATOMIC_COUNTER_GCC_41_H
#define PAGMO_ATOMIC_COUNTER_GCC_41_H

#include "base_atomic_counter.h"

namespace PaGMO
{
	template <class IntType>
	class atomic_counter_gcc_41: public base_atomic_counter<IntType, atomic_counter_gcc_41<IntType> >
	{
			typedef base_atomic_counter<IntType, atomic_counter_gcc_41<IntType> > ancestor;
		public:
			atomic_counter_gcc_41():ancestor::base_atomic_counter() {}
			template <class IntType2>
			atomic_counter_gcc_41 &operator+=(const IntType2 &n) {
				__sync_add_and_fetch(&(this->m_value),static_cast<IntType>(n));
				return *this;
			}
			template <class IntType2>
			atomic_counter_gcc_41 &operator-=(const IntType2 &n) {
				__sync_sub_and_fetch(&(this->m_value),static_cast<IntType>(n));
				return *this;
			}
	};
}

#endif
