# PaGMO lib source files.
SET(PAGMO_LIB_SRC_LIST
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/base.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/archipelago.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/base_island.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/island.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/population.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/sea.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/de.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/jde.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/mde_pbx.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/de_1220.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/pso.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/pso_generational.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/pso_generational_racing.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/bee_colony.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/firefly.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/cmaes.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/inverover.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/nn_tsp.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/nsga2.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/moea_d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/sms_emoa.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/cs.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/ihs.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/sga.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/sga_gray.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/monte_carlo.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/cstrs_co_evolution.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/null.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/sa_corana.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/mbh.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/ms.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/vega.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/cstrs_self_adaptive.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/pade.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/cstrs_immune_system.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/cstrs_core.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/nspso.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/spea2.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/migration/base.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/migration/base_r_policy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/migration/base_s_policy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/migration/best_s_policy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/migration/random_s_policy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/migration/best_kill_s_policy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/migration/hv_greedy_s_policy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/migration/hv_best_s_policy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/migration/fair_r_policy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/migration/worst_r_policy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/migration/random_r_policy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/migration/hv_greedy_r_policy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/migration/hv_fair_r_policy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/base.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/base_tsp.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/base_stochastic.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/base_dtlz.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/base_unc_mo.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/griewank.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/levy5.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/lennard_jones.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/lavor_maculan.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/ackley.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/branin.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/dtlz.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/golomb_ruler.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/identity.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/inventory.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/luksan_vlcek_1.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/luksan_vlcek_2.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/luksan_vlcek_3.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/himmelblau.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/bukin_f6.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/rastrigin.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/rosenbrock.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/schwefel.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/snopt_toyprob.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/string_match.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/tsp.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/tsp_cs.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/tsp_vrplc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/michalewicz.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/dejong.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/sch.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/fon.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/pol.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/kur.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/zdt.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/cec2006.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/cec2009.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/cec2013.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/con2mo.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/pressure_vessel.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/welded_beam.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/tens_comp_string.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/cstrs_co_evolution.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/death_penalty.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/cstrs_self_adaptive.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/antibodies_problem.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/shifted.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/scaled.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/rotated.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/normalized.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/decompose.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/noisy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/robust.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/con2uncon.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/barabasi_albert.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/clustered_ba.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/ageing_clustered_ba.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/base.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/custom.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/erdos_renyi.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/fully_connected.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/one_way_ring.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/pan.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/rim.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/ring.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/hypercube.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/unconnected.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/watts_strogatz.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/rng.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/util/hypervolume.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/util/hv_algorithm/base.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/util/hv_algorithm/hv2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/util/hv_algorithm/hv3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/util/hv_algorithm/wfg.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/util/hv_algorithm/bf_approx.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/util/hv_algorithm/bf_fpras.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/util/hv_algorithm/hoy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/util/hv_algorithm/hv4d_cpp_original/hv.c
	${CMAKE_CURRENT_SOURCE_DIR}/util/hv_algorithm/hv4d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/util/hv_algorithm/fpl_cpp_original/hv.c
	${CMAKE_CURRENT_SOURCE_DIR}/util/hv_algorithm/fpl.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/util/racing.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/util/discrepancy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/util/neighbourhood.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/util/race_pop.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/util/race_algo.cpp
)

# Additional files for the GTOP problems and keplerian toolbox.
IF(ENABLE_GTOP_DATABASE)
	SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST}
		${CMAKE_CURRENT_SOURCE_DIR}/problem/cassini_1.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/problem/cassini_2.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/problem/gtoc_1.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/problem/gtoc_2.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/problem/sagas.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/problem/rosetta.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/problem/messenger.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/problem/messenger_full.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/problem/tandem.cpp
	    ${CMAKE_CURRENT_SOURCE_DIR}/problem/tsp_ds.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/problem/laplace.cpp
#		${CMAKE_CURRENT_SOURCE_DIR}/problem/gtoc5_flyby.cpp
#		${CMAKE_CURRENT_SOURCE_DIR}/problem/gtoc5_launch.cpp
#		${CMAKE_CURRENT_SOURCE_DIR}/problem/gtoc5_rendezvous.cpp
#		${CMAKE_CURRENT_SOURCE_DIR}/problem/gtoc5_self_flyby.cpp
#		${CMAKE_CURRENT_SOURCE_DIR}/problem/earth_planet.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/problem/sample_return.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/problem/mga_1dsm_alpha.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/problem/mga_1dsm_tof.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/problem/mga_incipit.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/problem/mga_incipit_cstrs.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/problem/mga_part.cpp
#		${CMAKE_CURRENT_SOURCE_DIR}/problem/mga_target_event.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/mga.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/mga_dsm.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/misc4Tandem.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/Pl_Eph_An.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/PowSwingByInv.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/Lambert.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/Astro_Functions.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/time2distance.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/propagateKEP.cpp
	)
ENDIF(ENABLE_GTOP_DATABASE)

# Additional files for GSL support.
IF(ENABLE_GSL)
	SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST}
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/base_gsl.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/gsl_bfgs.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/gsl_bfgs2.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/gsl_derivative_free.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/gsl_fr.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/gsl_pr.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/gsl_gradient.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/gsl_nm.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/gsl_nm2.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/gsl_nm2rand.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/problem/spheres.cpp
#		${CMAKE_CURRENT_SOURCE_DIR}/problem/spheres_q.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/gsl_init.cpp
	)
ENDIF(ENABLE_GSL)

# Additional files for NLopt support.
IF(ENABLE_NLOPT)
	SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST}
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/base_nlopt.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/nlopt_cobyla.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/nlopt_bobyqa.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/nlopt_sbplx.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/nlopt_slsqp.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/nlopt_mma.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/nlopt_aug_lag.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/nlopt_aug_lag_eq.cpp
	)
ENDIF(ENABLE_NLOPT)

# Additional files for SNOPT support.
IF(ENABLE_SNOPT)
	SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST}
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/snopt.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/snopt_cpp_wrapper/snoptProblem_PAGMO.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/snopt_cpp_wrapper/snfilewrapper_PAGMO.c
	)
ENDIF(ENABLE_SNOPT)

# Additional files for IPOPT support.
IF(ENABLE_IPOPT)
	SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST}
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/ipopt.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/ipopt_cpp_wrapper/ipopt_problem.cpp
	)
ENDIF(ENABLE_IPOPT)

# Additional files for WORHP support.
IF(ENABLE_WORHP)
	SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST}
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/worhp.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/algorithm/worhp_cpp_wrapper/set_get_params.cpp
	)
ENDIF(ENABLE_WORHP)

# PaGMO MPI bits
IF(ENABLE_MPI)
	SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST}
		${CMAKE_CURRENT_SOURCE_DIR}/mpi_environment.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/mpi_island.cpp
	)
ENDIF(ENABLE_MPI)

# Create a pagmo_static library if main or tests or examples are requested.
IF(BUILD_MAIN OR ENABLE_TESTS OR BUILD_EXAMPLES)
	ADD_LIBRARY(pagmo_static STATIC ${PAGMO_LIB_SRC_LIST})
	SET_TARGET_PROPERTIES(pagmo_static PROPERTIES COMPILE_FLAGS "${STATIC_LIB_PAGMO_BUILD_FLAGS}" OUTPUT_NAME pagmo)
        INSTALL(TARGETS pagmo_static ARCHIVE DESTINATION ${LIB_INSTALL_PATH})
ENDIF(BUILD_MAIN OR ENABLE_TESTS OR BUILD_EXAMPLES)

# Create a pagmo shared library if PyGMO is requested.
IF(BUILD_PYGMO)
	ADD_LIBRARY(pagmo SHARED ${PAGMO_SHARED_LIB_SRC_LIST} ${PAGMO_LIB_SRC_LIST})
	# Request linking against the Python libs, in addition to our mandatory libs.
	TARGET_LINK_LIBRARIES(pagmo ${MANDATORY_LIBRARIES})
	# Add the definitions needed when building the PaGMO DLL.
	IF(ENABLE_GTOP_DATABASE)
		SET(DYNAMIC_LIB_PAGMO_BUILD_FLAGS "${DYNAMIC_LIB_PAGMO_BUILD_FLAGS} -DKEP_TOOL_DLL_EXPORT_API")
	ENDIF(ENABLE_GTOP_DATABASE)
	SET(DYNAMIC_LIB_PAGMO_BUILD_FLAGS "${DYNAMIC_LIB_PAGMO_BUILD_FLAGS} -DPAGMO_BUILD_PYGMO -DPAGMO_DLL_EXPORT_API")
	SET_TARGET_PROPERTIES(pagmo PROPERTIES COMPILE_FLAGS "${DYNAMIC_LIB_PAGMO_BUILD_FLAGS}")

	INSTALL(TARGETS pagmo
		RUNTIME DESTINATION ${LIB_INSTALL_PATH}
		LIBRARY DESTINATION ${LIB_INSTALL_PATH}
	)
ENDIF(BUILD_PYGMO)

IF(INSTALL_HEADERS)
	FILE(GLOB HEADERS_IN_EIGEN "${CMAKE_CURRENT_SOURCE_DIR}/Eigen/*")
	LIST(REMOVE_ITEM HEADERS_IN_EIGEN "${CMAKE_CURRENT_SOURCE_DIR}/Eigen/src")
	LIST(REMOVE_ITEM HEADERS_IN_EIGEN "${CMAKE_CURRENT_SOURCE_DIR}/Eigen/CMakeLists.txt")
	INSTALL(FILES ${HEADERS_IN_EIGEN} DESTINATION "include/pagmo/src/Eigen")
	INSTALL(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" DESTINATION "include/pagmo" FILES_MATCHING PATTERN "*.h")
	INSTALL(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../PyGMO" DESTINATION "include/pagmo" FILES_MATCHING PATTERN "*.h")
	INSTALL(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" DESTINATION "include/pagmo" FILES_MATCHING PATTERN "*.hpp")
	INSTALL(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../PyGMO" DESTINATION "include/pagmo" FILES_MATCHING PATTERN "*.hpp")
ENDIF(INSTALL_HEADERS)
