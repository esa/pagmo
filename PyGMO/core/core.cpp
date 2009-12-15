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

// 27/12/2008: Initial version by Francesco Biscani.

#include <boost/cstdint.hpp>
#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/manage_new_object.hpp>
#include <boost/python/module.hpp>
#include <boost/python/overloads.hpp>
#include <boost/utility.hpp>
#include <sstream>
#include <string>
#include <vector>

#include "../../src/GOclasses/algorithms/base.h"
#include "../../src/GOclasses/basic/archipelago.h"
#include "../../src/GOclasses/basic/topology/base_topology.h"
#include "../../src/GOclasses/basic/individual.h"
#include "../../src/GOclasses/basic/island.h"
#include "../../src/GOclasses/basic/population.h"
#include "../../src/GOclasses/problems/base.h"
#include "../../src/GOclasses/basic/migration/MigrationScheme.h"
#include "../../src/GOclasses/basic/migration/MigrationPolicy.h"
#include "../../src/GOclasses/basic/migration/MigrationSelectionPolicy.h"
#include "../../src/GOclasses/basic/migration/RandomMigrationSelectionPolicy.h"
#include "../../src/GOclasses/basic/migration/ChooseBestMigrationSelectionPolicy.h"
#include "../../src/GOclasses/basic/migration/MigrationReplacementPolicy.h"
#include "../../src/GOclasses/basic/migration/RandomMigrationReplacementPolicy.h"
#include "../../src/GOclasses/basic/migration/BestReplaceWorstIfBetterMigrationReplacementPolicy.h"
#include "../boost_python_container_conversions.h"
#include "../exceptions.h"
#include "../utils.h"

using namespace boost::python;
using namespace pagmo;

template <class T, class C>
static inline T *problem_getter(const C &c)
{
	return c.problem().clone();
}

template <class T, class C>
static inline T *algorithm_getter(const C &c)
{
	return c.algorithm().clone();
}

template <class T, class C>
static inline T *topology_getter(const C &c)
{
	return c.getTopology().clone();
}

/// \todo Is this really the correct way to do things??
template <class T, class C>
static inline T *selection_policy_getter(const C &c)
{
	return c.getMigrationSelectionPolicy().clone();
}

template <class T, class C>
static inline T *replacement_policy_getter(const C &c)
{
	return c.getMigrationReplacementPolicy().clone();
}

template <class T, class C>
static inline T *migration_scheme_getter(const C &c)
{
	return c.getMigrationScheme().clone();
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(island_evolve_overloads, evolve, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(archipelago_evolve_overloads, evolve, 0, 1)

// Instantiate the core module.
BOOST_PYTHON_MODULE(_core)
{
	// Translate exceptions for this module.
	translate_exceptions();

	// Enable automatic conversion to/from python list of double/vector<double>.
	to_tuple_mapping<std::vector<double> >();
	from_python_sequence<std::vector<double>,variable_capacity_policy>();

	// Enable automatic conversion to/from python list of double/vector<int>.
	to_tuple_mapping<std::vector<int> >();
	from_python_sequence<std::vector<int>,variable_capacity_policy>();

        // Enable automatic conversion to/from python list of double/vector<double>.
        to_tuple_mapping<std::vector<std::vector<double> > >();
        from_python_sequence<std::vector<std::vector<double> >,variable_capacity_policy>();

	// Expose individual class.
	class_<individual> class_ind("individual", "Individual class.", init<const problem::base &>());
	class_ind.def(init<const problem::base &, const std::vector<double> &, const std::vector<double> &>());
	class_ind.def(init<const problem::base &, const std::vector<double> &>());
	class_ind.def(init<const individual &>());
	class_ind.def("__copy__", &Py_copy_from_ctor<individual>);
	class_ind.def("__repr__", &Py_repr_from_stream<individual>);
	class_ind.add_property("fitness", &individual::getFitness, "Fitness.");
	class_ind.add_property("decision_vector", make_function(&individual::getDecisionVector, return_value_policy<copy_const_reference>()),
		"Decision vector.");
	class_ind.add_property("velocity", make_function(&individual::getVelocity, return_value_policy<copy_const_reference>()),
		"Velocity.");

	// Expose population class.
	typedef const individual &(population::*pop_get_const)(int) const;
	class_<population> class_pop("population", "Population class.", init<const problem::base &>());
	class_pop.def(init<const problem::base &, int>());
	class_pop.def(init<const problem::base &>());
	class_pop.def(init<const population &>());
	class_pop.def("__copy__", &Py_copy_from_ctor<population>);
	class_pop.def("__delitem__", &population::erase);
	class_pop.def("__getitem__", pop_get_const(&population::operator[]), return_value_policy<copy_const_reference>(), "Get a copy of individual.");
	class_pop.def("__len__", &population::size);
	class_pop.def("__setitem__", &population::setIndividual);
	class_pop.def("__repr__", &Py_repr_from_stream<population>);
	class_pop.add_property("problem", make_function(&problem_getter<problem::base,population>,return_value_policy<manage_new_object>()), "Problem.");
	class_pop.def("append", &population::push_back, "Append individual at the end of the population.");
	class_pop.def("insert", &population::insert, "Insert individual before index.");
	class_pop.def("mean", &population::evaluateMean, "Evaluate mean.");
	class_pop.def("std", &population::evaluateStd, "Evaluate std.");
	class_pop.def("best", &population::extractBestIndividual, return_value_policy<copy_const_reference>(), "Copy of best individual.");
	class_pop.def("worst", &population::extractWorstIndividual, return_value_policy<copy_const_reference>(), "Copy of worst individual.");

	// Expose island.
	class_<island> class_island("island", "Island.", init<const problem::base &, const algorithm::base &, int>());
	class_island.def(init<const problem::base &, const algorithm::base &>());
	class_island.def(init<const problem::base&, const algorithm::base&, int, const MigrationPolicy&>());
	class_island.def(init<const island &>());
	class_island.def("__copy__", &Py_copy_from_ctor<island>);
	class_island.def("__delitem__", &island::erase);
	class_island.def("__getitem__", &island::operator[], "Get a copy of individual.");
	class_island.def("__len__", &island::size);
	class_island.def("__setitem__", &island::set_individual);
	class_island.def("__repr__", &Py_repr_from_stream<island>);
	class_island.def("append", &island::push_back, "Append individual at the end of the island.");
	class_island.def("insert", &island::insert, "Insert individual after index.");
	class_island.add_property("problem", make_function(&problem_getter<problem::base,island>, return_value_policy<manage_new_object>()), "Problem.");
	class_island.add_property("algorithm", make_function(&algorithm_getter<algorithm::base,island>, return_value_policy<manage_new_object>()),
		&island::set_algorithm, "Algorithm.");
	class_island.add_property("selection_policy", make_function(&selection_policy_getter<MigrationSelectionPolicy, island>, return_value_policy<manage_new_object>()),
		&island::setMigrationSelectionPolicy, "The island's migration selection policy.");
	class_island.add_property("replacement_policy", make_function(&replacement_policy_getter<MigrationReplacementPolicy, island>, return_value_policy<manage_new_object>()),
		&island::setMigrationReplacementPolicy, "The island's migration selection policy.");
	class_island.add_property("population", &island::get_population, "Copy of population.");
	class_island.def("mean", &island::mean, "Evaluate mean.");
	class_island.def("std", &island::std, "Evaluate std.");
	class_island.def("best", &island::best, "Copy of best individual.");
	class_island.def("worst", &island::worst, "Copy of worst individual.");
	class_island.add_property("id", &island::id, "Identification number.");
	class_island.def("evolve", &island::evolve, island_evolve_overloads());
	class_island.def("evolve_t", &island::evolve_t, "Evolve for an amount of time.");
	class_island.def("join", &island::join, "Block until evolution has terminated.");
	class_island.add_property("busy", &island::busy, "True if island is evolving, false otherwise.");
	class_island.add_property("evo_time", &island::evo_time, "Total time spent evolving.");

	// Expose archipelago.
	class_<archipelago> class_arch("archipelago", "Archipelago", init<const problem::base &>());
	class_arch.def(init<const problem::base &, const algorithm::base &, int, int>());
	class_arch.def(init<const problem::base &, const MigrationScheme&>());
	class_arch.def(init<const problem::base &, const algorithm::base &, int, int, const Migration&>());
	class_arch.def(init<const archipelago &>());
	class_arch.def("__copy__", &Py_copy_from_ctor<archipelago>);
	class_arch.def("__getitem__", &archipelago::operator[], return_value_policy<copy_const_reference>());
	class_arch.def("__len__", &archipelago::size);
	class_arch.def("__setitem__", &archipelago::set_island);
	class_arch.def("__repr__", &Py_repr_from_stream<archipelago>);
	class_arch.add_property("migration_scheme", make_function(&migration_scheme_getter<MigrationScheme, archipelago>, return_value_policy<manage_new_object>()),
		&archipelago::setMigrationScheme, "The archipelago's migration scheme.");
	class_arch.add_property("topology", make_function(&topology_getter<base_topology, archipelago>, return_value_policy<manage_new_object>()),
		&archipelago::setTopology, "The archipelago's migration topology.");
	class_arch.def("append", &archipelago::push_back, "Append island.");
	class_arch.add_property("problem", make_function(&problem_getter<problem::base,archipelago>, return_value_policy<manage_new_object>()), "Problem.");
	class_arch.def("join", &archipelago::join, "Block until evolution on each island has terminated.");
	class_arch.add_property("busy", &archipelago::busy, "True if at least one island is evolving, false otherwise.");
	class_arch.def("evolve", &archipelago::evolve, archipelago_evolve_overloads());
	class_arch.def("evolve_t", &archipelago::evolve_t, "Evolve islands for an amount of time.");
	class_arch.def("best", &archipelago::best, "Copy of best individual.");
	class_arch.def("max_evo_time", &archipelago::getMaxEvoTime, "Maximum of total evolution times for all islands.");
	class_arch.def("total_evo_time", &archipelago::getTotalEvoTime, "Sum of total evolution times for all islands.");
	
	
	// Expose Migration
	class_<Migration> class_M("migration", "The migration parameters for archipelago.", init<>());
	class_M.def(init<const MigrationScheme&, const MigrationPolicy&>());
	class_M.def("__copy__", &Py_copy_from_ctor<Migration>);
	class_M.def("__repr__", &Py_repr_from_stream<Migration>);

	// Expose MigrationScheme
	class_<MigrationScheme> class_MS("migration_scheme", "The migration scheme.", init<int, int, optional<boost::uint32_t> >());
	class_MS.def(init<int, int, const base_topology&, optional<boost::uint32_t> >());
	class_MS.def("__copy__", &Py_copy_from_ctor<MigrationScheme>);
	class_MS.def("__repr__", &Py_repr_from_stream<MigrationScheme>);
	
	// Expose MigrationPolicy
	class_<MigrationPolicy> class_MP("migration_policy", "The island's migration policy.", init<>());
	class_MP.def(init<const double>());
	class_MP.def(init<const double, const MigrationSelectionPolicy&, const MigrationReplacementPolicy&>());
	class_MP.def("__copy__", &Py_copy_from_ctor<MigrationPolicy>);
	class_MP.def("__repr__", &Py_repr_from_stream<MigrationPolicy>);
	
	
	// Expose MigrationSelectionPolicy
	class_<MigrationSelectionPolicy, boost::noncopyable> class_MSP("__migration_selection_policy", "A migration selection policy.", no_init);
	
	// Expose RandomMigrationSelectionPolicy
	class_<RandomMigrationSelectionPolicy, bases<MigrationSelectionPolicy> > class_RMSP("random_selection_policy", "A random migration selection policy.", init<optional<const boost::uint32_t> >());	
	/*
	 * !!!
	 * Here and below the order of declaration of constructors is crucial !!!
	 * If the int one is given first, RandomMigrationSelectionPolicy(1) will be interpreted as RandomMigrationSelectionPolicy(1.0) !!!
	 * Python overloading SUCKS.
	 * !!!
	 */	
	class_RMSP.def(init<const double&, optional<const boost::uint32_t> >());
	class_RMSP.def(init<const int&, optional<const boost::uint32_t> >());
	class_RMSP.def("__copy__", &Py_copy_from_ctor<RandomMigrationSelectionPolicy>);
	class_RMSP.def("__repr__", &Py_repr_from_stream<RandomMigrationSelectionPolicy>);
	// Expose ChooseBestMigrationSelectionPolicy
	class_<ChooseBestMigrationSelectionPolicy, bases<MigrationSelectionPolicy> > class_CBMSP("choose_best_selection_policy", "A choose best migration selection policy.", init<>());	
	class_CBMSP.def(init<const double&>());
	class_CBMSP.def(init<const int&>());	
	class_CBMSP.def("__copy__", &Py_copy_from_ctor<ChooseBestMigrationSelectionPolicy>);
	class_CBMSP.def("__repr__", &Py_repr_from_stream<ChooseBestMigrationSelectionPolicy>);
		
	// Expose MigrationReplacementPolicy
	class_<MigrationReplacementPolicy, boost::noncopyable> class_MRP("__migration_replacement_policy", "A migration replacement policy.", no_init);
	
	// Expose RandomMigrationReplacementPolicy
	class_<RandomMigrationReplacementPolicy, bases<MigrationReplacementPolicy> > class_RMRP("random_replacement_policy", "A random migration replacement policy.", init<optional<const boost::uint32_t> >());
	class_RMRP.def(init<const double&, optional<const boost::uint32_t> >());
	class_RMRP.def(init<const int&, optional<const boost::uint32_t> >());	
	class_RMRP.def("__copy__", &Py_copy_from_ctor<RandomMigrationReplacementPolicy>);
	class_RMRP.def("__repr__", &Py_repr_from_stream<RandomMigrationReplacementPolicy>);

	class_<BestReplaceWorstIfBetterMigrationReplacementPolicy, bases<MigrationReplacementPolicy> > class_BRWIBMRP("best_rep_worst_replacement_policy", "A migration replacement policy where best incoming replace worst present if they are better.", init<>());
	class_BRWIBMRP.def(init<const double&>());
	class_BRWIBMRP.def(init<const int&>());
	class_BRWIBMRP.def("__copy__", &Py_copy_from_ctor<BestReplaceWorstIfBetterMigrationReplacementPolicy>);
	class_BRWIBMRP.def("__repr__", &Py_repr_from_stream<BestReplaceWorstIfBetterMigrationReplacementPolicy>);	
}
