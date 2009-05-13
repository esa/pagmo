# PaGMO lib source files.
SET(PAGMO_LIB_SRC_LIST
	${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/Astro_Functions.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/Lambert.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/mga.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/mga_dsm.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/misc4Tandem.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/Pl_Eph_An.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/PowSwingByInv.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/propagateKEP.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/time2distance.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Functions/objfuns/classicobjfuns.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Functions/objfuns/trajobjfuns.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/algorithms/ASA.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/algorithms/DE.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/algorithms/MPSO.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/algorithms/PSO.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/algorithms/SGA.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/algorithms/go_algorithm.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/algorithms/ihs_algorithm.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/algorithms/nm_algorithm.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/individual.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/archipelago.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/ba_topology.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/base_topology.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/fully_connected_topology.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/graph_topology.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/one_way_ring_topology.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/chain_topology.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/broadcast_topology.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/ring_topology.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/ring12_topology.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/ring123_topology.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/torus_topology.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/cartwheel_topology.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/lattice_topology.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/hypercube_topology.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/island.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/population.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/MigrationScheme.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/MigrationSelectionPolicy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/RandomMigrationSelectionPolicy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/ChooseBestMigrationSelectionPolicy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/MigrationReplacementPolicy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/RandomMigrationReplacementPolicy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/BestReplaceWorstIfBetterMigrationReplacementPolicy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/problems/ClassicProblems.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/problems/GOproblem.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/problems/TrajectoryProblems.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/problems/twodee_problem.cpp
)

IF(UNIX)
	SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST}
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/pthread/exceptions.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/pthread/once.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/pthread/thread.cpp
	)
ELSE(UNIX)
	SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST}
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/win32/exceptions.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/win32/thread.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/win32/tss_dll.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/win32/tss_pe.cpp
	)
ENDIF(UNIX)

# Create a pagmo_static library if main is requested.
IF(BUILD_MAIN)
	SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST}
		${CMAKE_CURRENT_SOURCE_DIR}/SolversThreads/SolversThreads.cpp
	)
	ADD_LIBRARY(pagmo_static STATIC ${PAGMO_LIB_SRC_LIST})
	SET_TARGET_PROPERTIES(pagmo_static PROPERTIES COMPILE_FLAGS "-DBOOST_THREAD_BUILD_LIB")
ENDIF(BUILD_MAIN)

# Create a pagmo shared library if PyGMO is requested.
IF(BUILD_PYGMO)
	# We need to compile the source of the Boost.Python libraries.
	SET(PAGMO_SHARED_LIB_SRC_LIST
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/converter/arg_to_python_base.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/converter/builtin_converters.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/converter/from_python.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/converter/registry.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/converter/type_id.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/object/class.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/object/enum.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/object/function.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/object/function_doc_signature.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/object/inheritance.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/object/iterator.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/object/life_support.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/object/pickle_support.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/object/stl_iterator.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/dict.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/errors.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/exec.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/import.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/list.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/long.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/module.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/numeric.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/object_operators.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/object_protocol.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/slice.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/str.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/tuple.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_python/wrapper.cpp
	)
	ADD_LIBRARY(pagmo SHARED ${PAGMO_SHARED_LIB_SRC_LIST} ${PAGMO_LIB_SRC_LIST})
	# Request linking against the Python libs, in addition to our mandatory libs.
	TARGET_LINK_LIBRARIES(pagmo ${MANDATORY_LIBRARIES} ${PYTHON_LIBRARY})
	# A couple of extra definitions. The first is required when building Boost.Python,
	# the second is to signal the preprocessor that we are building PyGMO (used, e.g., in exceptions.h).
	IF(UNIX)
		SET_TARGET_PROPERTIES(pagmo PROPERTIES COMPILE_FLAGS "-DBOOST_PYTHON_SOURCE -DBOOST_THREAD_BUILD_DLL -DPAGMO_BUILD_PYGMO -DPAGMO_DLL_EXPORT_API")
	ELSE(UNIX)
		SET_TARGET_PROPERTIES(pagmo PROPERTIES COMPILE_FLAGS "-DBOOST_PYTHON_SOURCE -DBOOST_THREAD_BUILD_DLL -DPAGMO_BUILD_PYGMO -DPAGMO_DLL_EXPORT_API -DBOOST_HAS_DECLSPEC")
	ENDIF(UNIX)

	INSTALL(TARGETS pagmo
		RUNTIME DESTINATION ${LIB_INSTALL_PATH}
		LIBRARY DESTINATION ${LIB_INSTALL_PATH}
	)
ENDIF(BUILD_PYGMO)
