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

#ifndef PAGMO_ALGORITHMS_H
#define PAGMO_ALGORITHMS_H

// Header including all algorithms implemented in PaGMO.

// Heuristics
#include "algorithm/base.h"
#include "algorithm/cs.h"
#include "algorithm/de.h"
#include "algorithm/de_1220.h"
#include "algorithm/sea.h"
#include "algorithm/jde.h"
#include "algorithm/ihs.h"
#include "algorithm/mde_pbx.h"
#include "algorithm/monte_carlo.h"
#include "algorithm/null.h"
#include "algorithm/pso.h"
#include "algorithm/pso_generational.h"
#include "algorithm/pso_generational_racing.h"
#include "algorithm/sa_corana.h"
#include "algorithm/sga.h"
#include "algorithm/sga_gray.h"
#include "algorithm/nsga2.h"
#include "algorithm/bee_colony.h"
#include "algorithm/firefly.h"
#include "algorithm/cmaes.h"
#include "algorithm/nsga2.h"
#include "algorithm/vega.h"
#include "algorithm/cstrs_co_evolution.h"
#include "algorithm/pade.h"
#include "algorithm/moea_d.h"
#include "algorithm/cstrs_immune_system.h"
#include "algorithm/cstrs_core.h"
#include "algorithm/sms_emoa.h"
#include "algorithm/nspso.h"
#include "algorithm/spea2.h"
#include "algorithm/inverover.h"
#include "algorithm/nn_tsp.h"

// Hyper-heuristics
#include "algorithm/mbh.h"
#include "algorithm/ms.h"
#include "algorithm/cstrs_self_adaptive.h"

// SNOPT algorithm.
#ifdef PAGMO_ENABLE_SNOPT
	#include "algorithm/snopt.h"
#endif

// IPOPT algorithm.
#ifdef PAGMO_ENABLE_IPOPT
	#include "algorithm/ipopt.h"
#endif

// GSL algorithms.
#ifdef PAGMO_ENABLE_GSL
	#include "algorithm/base_gsl.h"
	#include "algorithm/gsl_bfgs.h"
	#include "algorithm/gsl_bfgs2.h"
	#include "algorithm/gsl_fr.h"
	#include "algorithm/gsl_nm.h"
	#include "algorithm/gsl_nm2.h"
	#include "algorithm/gsl_nm2rand.h"
	#include "algorithm/gsl_pr.h"
#endif

// NLopt algorithms.
#ifdef PAGMO_ENABLE_NLOPT
	#include "algorithm/nlopt_bobyqa.h"
	#include "algorithm/nlopt_cobyla.h"
	#include "algorithm/nlopt_sbplx.h"
	#include "algorithm/nlopt_slsqp.h"
	#include "algorithm/nlopt_mma.h"
	#include "algorithm/nlopt_aug_lag.h"
	#include "algorithm/nlopt_aug_lag_eq.h"
#endif

#ifdef PAGMO_ENABLE_WORHP
	#include "algorithm/worhp.h"
#endif

#endif
