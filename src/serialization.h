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

#ifndef PAGMO_SERIALIZATION_H
#define PAGMO_SERIALIZATION_H

#include "adj_list_serialize_bug_fix.hpp"

// Headers needed for serialization purposes.
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
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
#include <limits>

// Serialization of circular buffer, unordered map.
// TODO: serialize the functors.. allocator, Hash, Pred, etc.
#include <boost/circular_buffer.hpp>
#include <boost/unordered_map.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <string>
#include <utility>
#include <vector>

#include "Eigen/Dense"

namespace pagmo {

/// Custom save function for the serialization of vector of doubles that handle also inf and NaN.
template <class Archive>
void custom_vector_double_save(Archive &ar, const std::vector<double> &v, const unsigned int)
{
	const std::vector<double>::size_type size = v.size();
	// Save size.
	ar << size;
	// Save elements.
	std::string tmp;
	for (std::vector<double>::size_type i = 0; i < size; ++i) {
		if (boost::math::isnan(v[i])) {
			tmp = "nan";
		} else if (boost::math::isinf(v[i])) {
			if (v[i] > 0) {
				tmp = "inf";
			} else {
				tmp = "-inf";
			}
		} else {
			tmp = boost::lexical_cast<std::string>(v[i]);
		}
		ar << tmp;
	}
}

/// Custom load function for the serialization of vector of doubles that handle also inf and NaN.
template <class Archive>
void custom_vector_double_load(Archive &ar, std::vector<double> &v, const unsigned int)
{
	std::vector<double>::size_type size = 0;
	// Load size.
	ar >> size;
	v.resize(size);
	std::string tmp;
	for (std::vector<double>::size_type i = 0; i < size; ++i) {
		ar >> tmp;
		if (tmp == "nan") {
			v[i] = std::numeric_limits<double>::quiet_NaN();
		} else if (tmp == "inf") {
			v[i] = std::numeric_limits<double>::infinity();
		} else if (tmp == "-inf") {
			v[i] = -std::numeric_limits<double>::infinity();
		} else {
			v[i] = boost::lexical_cast<double>(tmp);
		}
	}
}

}

namespace boost { namespace serialization {

template <class Archive>
void save(Archive &ar, const boost::circular_buffer<std::vector<double> > &cb, unsigned int version)
{
	// Let's first save capacity and size.
	typedef typename boost::circular_buffer<std::vector<double> >::capacity_type capacity_type;
	typedef typename boost::circular_buffer<std::vector<double> >::size_type size_type;
	capacity_type capacity = cb.capacity();
	ar << capacity;
	size_type size = cb.size();
	ar << size;
	// Save the content.
	for (size_type i = 0; i < size; ++i) {
		pagmo::custom_vector_double_save(ar,cb[i],version);
	}
}

template <class Archive>
void load(Archive &ar, boost::circular_buffer<std::vector<double> > &cb, unsigned int version)
{
	typedef typename boost::circular_buffer<std::vector<double> >::capacity_type capacity_type;
	typedef typename boost::circular_buffer<std::vector<double> >::size_type size_type;
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
		pagmo::custom_vector_double_load(ar,cb[i],version);
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

template <class Archive>
void save(Archive &ar, const Eigen::MatrixXd &cb, unsigned int)
{
	// Let's first save the dimension
	std::size_t nrows = cb.rows(), ncols = cb.cols();
	ar << nrows;
	ar << ncols;
	//And then the numbers
	for (unsigned int i = 0; i < (unsigned int)cb.rows(); ++i) {
		for (unsigned int j = 0; j < (unsigned int)cb.cols(); ++j) {
			ar << cb(i,j);
		}
	}
}

template <class Archive>
void load(Archive &ar, Eigen::MatrixXd &cb, unsigned int)
{
	std::size_t rows, cols;
	// Let's first restore the dimension
	ar >> rows;
	ar >> cols;
	cb.resize(rows,cols);
	//And then the numbers
	for (unsigned int i = 0; i < (unsigned int)cb.rows(); ++i) {
		for (std::size_t j = 0; j < (unsigned int)cb.cols(); ++j) {
			ar >> cb(i,j);
		}
	}
}

template <class Archive>
void serialize(Archive &ar, Eigen::MatrixXd &cb, const unsigned int version)
{
	split_free(ar, cb, version);
}

template <class Archive>
void save(Archive &ar, const Eigen::VectorXd &cb, unsigned int)
{
	std::size_t nrows = cb.rows();
	// Let's first save the dimension
	ar << nrows;
	//And then the numbers
	for (unsigned int i = 0; i < (unsigned int)cb.rows(); ++i) {
		ar << cb(i);
	}
}

template <class Archive>
void load(Archive &ar, Eigen::VectorXd &cb, unsigned int)
{
	std::size_t rows;
	// Let's first restore the dimension
	ar >> rows;
	cb.resize(rows);
	//And then the numbers
	for (unsigned int i = 0; i < (unsigned int)cb.rows(); ++i) {
		ar >> cb(i);
	}
}

template <class Archive>
void serialize(Archive &ar, Eigen::VectorXd &cb, const unsigned int version)
{
	split_free(ar, cb, version);
}

template <class Archive>
void serialize(Archive &ar, boost::vecS &cb, const unsigned int version) 
{ 
	(void)(ar);
	(void)(cb);
	(void)(version);
}


// ---------- boost::normal_distribution<double> ----------
// Serialization of boost::normal_distribution<double> exploits the fact that the
// state of Boost distributions can be sent/received to/from standard streams.
template <class Archive>
void save(Archive &ar, const boost::normal_distribution<double> &cb, const unsigned int)
{ 
	std::stringstream ss;
	ss << cb;
	std::string tmp(ss.str());
	ar << tmp;
}

template <class Archive>
void load(Archive &ar, boost::normal_distribution<double> &cb, const unsigned int)
{
	std::string tmp;
	ar >> tmp;
	std::stringstream ss(tmp);
	ss >> cb;
}

template <class Archive>
void serialize(Archive &ar, boost::normal_distribution<double> &cb, const unsigned int version)
{
	split_free(ar, cb, version);
}

// ---------- boost::uniform_real_distribution -----------
// Serialization of boost::uniform_real_distribution exploits the fact that the
// state of Boost distributions can be sent/received to/from standard streams.
template <class Archive>
void save(Archive &ar, const boost::random::uniform_real_distribution<double> &cb, const unsigned int)
{ 
	std::stringstream ss;
	ss << cb;
	std::string tmp(ss.str());
	ar << tmp;
}

template <class Archive>
void load(Archive &ar, boost::random::uniform_real_distribution<double> &cb, const unsigned int)
{
	std::string tmp;
	ar >> tmp;
	std::stringstream ss(tmp);
	ss >> cb;
}

template <class Archive>
void serialize(Archive &ar, boost::random::uniform_real_distribution<double> &cb, const unsigned int version)
{
	split_free(ar, cb, version);
}

}}

#endif
