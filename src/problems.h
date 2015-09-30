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

#ifndef PAGMO_PROBLEMS_H
#define PAGMO_PROBLEMS_H

// Header including all problems implemented in PaGMO.

#include "problem/ackley.h"
#include "problem/base.h"
#include "problem/base_stochastic.h"
#include "problem/base_dtlz.h"
#include "problem/base_unc_mo.h"
#include "problem/branin.h"
#include "problem/bukin_f6.h"
#include "problem/cec2006.h"
#include "problem/cec2009.h"
#include "problem/cec2013.h"
#include "problem/con2mo.h"
#include "problem/dejong.h"
#include "problem/dtlz.h"
#include "problem/golomb_ruler.h"
#include "problem/griewank.h"
#include "problem/himmelblau.h"
#include "problem/identity.h"
#include "problem/inventory.h"
#include "problem/luksan_vlcek_1.h"
#include "problem/luksan_vlcek_2.h"
#include "problem/luksan_vlcek_3.h"
#include "problem/lennard_jones.h"
#include "problem/lavor_maculan.h"
#include "problem/levy5.h"
#include "problem/michalewicz.h"
#include "problem/rastrigin.h"
#include "problem/rosenbrock.h"
#include "problem/schwefel.h"
#include "problem/sch.h"
#include "problem/snopt_toyprob.h"
#include "problem/string_match.h"
#include "problem/tsp.h"
#include "problem/tsp_vrplc.h"
#include "problem/tsp_cs.h"
#include "problem/fon.h"
#include "problem/pol.h"
#include "problem/kur.h"
#include "problem/zdt.h"
#include "problem/pressure_vessel.h"
#include "problem/tens_comp_string.h"
#include "problem/welded_beam.h"
#include "problem/death_penalty.h"
#include "problem/cstrs_co_evolution.h"
#include "problem/cstrs_self_adaptive.h"
#include "problem/antibodies_problem.h"
#include "problem/shifted.h"
#include "problem/scaled.h"
#include "problem/rotated.h"
#include "problem/normalized.h"
#include "problem/decompose.h"
#include "problem/noisy.h"
#include "problem/robust.h"
#include "problem/con2uncon.h"

// GSL problems.
#ifdef PAGMO_ENABLE_GSL
	#include "problem/spheres.h"
//	#include "problem/spheres_q.h"
#endif

#ifdef PAGMO_ENABLE_KEP_TOOLBOX
        #include "problem/cassini_1.h"
        #include "problem/cassini_2.h"
        #include "problem/gtoc_1.h"
        #include "problem/gtoc_2.h"
        #include "problem/sagas.h"
        #include "problem/rosetta.h"
        #include "problem/messenger.h"
        #include "problem/messenger_full.h"
        #include "problem/tandem.h"
        #include "problem/laplace.h"
//        #include "problem/sample_return.h"
//        #include "problem/gtoc5_flyby.h"
//        #include "problem/gtoc5_launch.h"
//        #include "problem/gtoc5_rendezvous.h"
//        #include "problem/gtoc5_self_flyby.h"
//        #include "problem/earth_planet.h"
        #include "problem/mga_1dsm_alpha.h"
        #include "problem/mga_1dsm_tof.h"
        #include "problem/mga_incipit.h"
	#include "problem/mga_incipit_cstrs.h"
        #include "problem/mga_part.h"
	#include "problem/mga_target_event.h"
        #include "problem/tsp_ds.h"
#endif

#endif
