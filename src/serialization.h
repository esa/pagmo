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

#ifndef PAGMO_SERIALIZATION_H
#define PAGMO_SERIALIZATION_H

// Headers needed for serialization purposes.
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

// Serialization of circular buffer, unordered map.
// TODO: serialize the functors.. allocator, Hash, Pred, etc.
#include <boost/circular_buffer.hpp>
#include <boost/unordered_map.hpp>
#include <utility>

namespace boost { namespace serialization {

template <class Archive, class T>
void save(Archive &ar, const boost::circular_buffer<T> &cb, unsigned int)
{
	// Let's first save capacity and size.
	typedef typename boost::circular_buffer<T>::capacity_type capacity_type;
	typedef typename boost::circular_buffer<T>::size_type size_type;
	capacity_type capacity = cb.capacity();
	ar << capacity;
	size_type size = cb.size();
	ar << size;
	// Save the content.
	for (size_type i = 0; i < size; ++i) {
		ar << cb[i];
	}
}

template <class Archive, class T>
void load(Archive &ar, boost::circular_buffer<T> &cb, unsigned int)
{
	typedef typename boost::circular_buffer<T>::capacity_type capacity_type;
	typedef typename boost::circular_buffer<T>::size_type size_type;
	// Restore capacity.
	capacity_type capacity;
	ar >> capacity;
	cb.set_capacity(capacity);
	// Restore size.
	size_type size;
	ar >> size;
	cb.resize(size);
	// Restore elements.
	for (size_type i = 0; i < size; ++i) {
		ar >> cb[i];
	}
}

template <class Archive, class T>
void serialize(Archive &ar, boost::circular_buffer<T> &cb, const unsigned int version)
{
	split_free(ar, cb, version);
}

template <class Archive, class Key, class Mapped, class Hash, class Pred, class Alloc>
void save(Archive &ar, const boost::unordered_map<Key,Mapped,Hash,Pred,Alloc> &um, unsigned int)
{
	typedef typename boost::unordered_map<Key,Mapped,Hash,Pred,Alloc>::size_type size_type;
	typedef typename boost::unordered_map<Key,Mapped,Hash,Pred,Alloc>::const_iterator const_iterator;
	// Save size.
	const size_type size = um.size();
	ar << size;
	for (const_iterator it = um.begin(); it != um.end(); ++it) {
		ar << *it;
	}
}

template <class Archive, class Key, class Mapped, class Hash, class Pred, class Alloc>
void load(Archive &ar, boost::unordered_map<Key,Mapped,Hash,Pred,Alloc> &um, unsigned int)
{
	typedef typename boost::unordered_map<Key,Mapped,Hash,Pred,Alloc>::size_type size_type;
	// Empty the map.
	um.clear();
	// Recover size.
	size_type size;
	ar >> size;
	for (size_type i = 0; i < size; ++i) {
		std::pair<Key, Mapped> value;
		ar >> value;
		um[value.first] = value.second;
	}
}

template <class Archive, class Key, class Mapped, class Hash, class Pred, class Alloc>
void serialize(Archive &ar, boost::unordered_map<Key,Mapped,Hash,Pred,Alloc> &um, unsigned int version)
{
	split_free(ar, um, version);
}

}}

#endif
