/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
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

#include <Python.h>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/utility.hpp>

#include "../../src/util/hypervolume.h"
#include "../../src/util/hv_algorithm/base.h"
#include "../../src/util/hv_algorithm/lebmeasure.h"
#include "../../src/util/hv_algorithm/optimal2d.h"
#include "../../src/util/hv_algorithm/optimal3d.h"
#include "../utils.h"

using namespace boost::python;
using namespace pagmo;

template <class HVAlgorithm>
static inline class_<HVAlgorithm,bases<util::hv_algorithm::base> > algorithm_wrapper(const char *name, const char *descr)
{

	class_<HVAlgorithm,bases<util::hv_algorithm::base> > retval(name,descr,init<const HVAlgorithm &>());
	retval.def(init<>());
	retval.def("compute", &HVAlgorithm::compute);
	return retval;
}

void expose_hv_algorithm() {
	class_<util::hv_algorithm::base,boost::noncopyable>("_base",no_init)
		.def("compute", &util::hv_algorithm::base::compute);
	algorithm_wrapper<util::hv_algorithm::lebmeasure>("lebmeasure","LebMeasure algorithm.");
	algorithm_wrapper<util::hv_algorithm::optimal2d>("optimal2d","Optimal2D algorithm.");
	algorithm_wrapper<util::hv_algorithm::optimal3d>("optimal3d","Optimal3D algorithm.");
}

void expose_hypervolume() {
	class_<util::hypervolume>("hypervolume","Hypervolume class.", init<const std::vector<std::vector<double> > &>())
		.def(init<boost::shared_ptr<population> >())
		.def("compute", &util::hypervolume::compute);
}

BOOST_PYTHON_MODULE(_util) {
	common_module_init();

	expose_hypervolume();

	scope current;
	std::string submoduleName(extract<const char*>(current.attr("__name__")));
	submoduleName.append(".hv_algorithm");

	object submodule(borrowed(PyImport_AddModule(submoduleName.c_str())));
	current.attr("hv_algorithm") = submodule;
	scope submoduleScope = submodule;
	expose_hv_algorithm();

}
