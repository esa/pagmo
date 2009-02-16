/* Copyright 2006-2008 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org/libs/flyweight for library home page.
 */

#ifndef BOOST_FLYWEIGHT_REFCOUNTED_HPP
#define BOOST_FLYWEIGHT_REFCOUNTED_HPP

#if defined(_MSC_VER)&&(_MSC_VER>=1200)
#pragma once
#endif

#include <boost/config.hpp> /* keep it first to prevent nasty warns in MSVC */
#include <algorithm>
#include <boost/assert.hpp>
#include <boost/detail/atomic_count.hpp>
#include <boost/detail/workaround.hpp>
#include <boost/flyweight/refcounted_fwd.hpp>
#include <boost/flyweight/tracking_tag.hpp>
#include <boost/utility/swap.hpp>

/* Refcounting tracking policy: values have an embedded ref count,
 * when this goes down to zero the element is erased from the
 * factory.
 */

namespace boost{

namespace flyweights{

namespace detail{

template<typename Value,typename Key>
class refcounted_value
{
public:
  explicit refcounted_value(const Value& x_):
    x(x_),ref(0)
  {}
  
  refcounted_value(const refcounted_value& r):
    x(r.x),ref(0)
  {}

  ~refcounted_value()
  {
    /* count()!=0 most likely indicates that the flyweight factory
     * has been destructed before some of the flyweight objects using
     * it. Check for static initialization order problems with this
     * flyweight type.
     */

    BOOST_ASSERT(count()==0);
  }

  refcounted_value& operator=(const refcounted_value& r)
  {
    x=r.x;
    return *this;
  }
  
  operator const Value&()const{return x;}
  operator const Key&()const{return x;}
    
#if !defined(BOOST_NO_MEMBER_TEMPLATE_FRIENDS)
private:
  template<typename,typename> friend class refcounted_handle;
#endif

  long count()const{return ref;}
  void add_ref()const{++ref;}
  bool release()const{return (--ref==0);}

private:
  Value                               x;
  mutable boost::detail::atomic_count ref;
};

template<typename Handle,typename TrackingHelper>
class refcounted_handle
{
public:
  explicit refcounted_handle(const Handle& h_):h(h_)
  {
    TrackingHelper::entry(*this).add_ref();
  }
  
  refcounted_handle(const refcounted_handle& x):h(x.h)
  {
    TrackingHelper::entry(*this).add_ref();
  }

  refcounted_handle& operator=(refcounted_handle x)
  {
    swap(*this,x);
    return *this;
  }

  ~refcounted_handle()
  {
    if(TrackingHelper::entry(*this).release()){
      TrackingHelper::erase(*this,check_erase);
    }
  }

  operator const Handle&()const{return h;}

  friend void swap(refcounted_handle& x, refcounted_handle& y)
  {
    boost::swap(x.h,y.h);
  }

private:
  static bool check_erase(const refcounted_handle& x)
  {
    return TrackingHelper::entry(x).count()==0;
  }

  Handle h;
};

} /* namespace flyweights::detail */

struct refcounted:tracking_marker
{
  struct entry_type
  {
    template<typename Value,typename Key>
    struct apply
    {
      typedef detail::refcounted_value<Value,Key> type;
    };
  };

  struct handle_type
  {
    template<typename Handle,typename TrackingHelper>
    struct apply
    {
      typedef detail::refcounted_handle<Handle,TrackingHelper> type;
    };
  };
};

} /* namespace flyweights */

} /* namespace boost */

#endif
