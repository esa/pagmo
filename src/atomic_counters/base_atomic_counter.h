#ifndef PAGMO_BASE_ATOMIC_COUNTER_H
#define PAGMO_BASE_ATOMIC_COUNTER_H

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace PaGMO
{
	template <class IntType, class Derived>
	class base_atomic_counter
	{
		public:
			base_atomic_counter():m_value(0) {}
			template <class IntType2>
			base_atomic_counter(const IntType2 &n):m_value(static_cast<IntType2>(n)) {}
			operator IntType() const {
				return m_value;
			}
		protected:
			IntType m_value;
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
