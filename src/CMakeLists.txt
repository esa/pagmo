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
	${CMAKE_CURRENT_SOURCE_DIR}/Functions/objfuns/trajobjfuns.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/algorithms/ASA.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/algorithms/DE.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/algorithms/MPSO.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/algorithms/PSO.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/algorithms/SGA.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/individual.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/population.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/problems/TrajectoryProblems.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SolversThreads/SolversThreads.cpp
)

# Boost thread setup.
ADD_DEFINITIONS(-DBOOST_THREAD_BUILD_LIB)
SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST} boost_thread/tss_null.cpp)
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
	ADD_LIBRARY(pagmo_static STATIC ${PAGMO_LIB_SRC_LIST})
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
	# This definition is required when compiling Boost.Python.
	ADD_DEFINITIONS(-DBOOST_PYTHON_SOURCE)
	ADD_LIBRARY(pagmo SHARED ${PAGMO_SHARED_LIB_SRC_LIST} ${PAGMO_LIB_SRC_LIST})
	# Request linking against the Python libs, in addition to our mandatory libs.
	TARGET_LINK_LIBRARIES(pagmo ${MANDATORY_LIBRARIES} ${PYTHON_LIBRARY})
	# Additional definition to signal the preprocessor that we are building PyGMO (used, e.g., in exceptions.h).
	SET_TARGET_PROPERTIES(pagmo PROPERTIES COMPILE_FLAGS "-DPAGMO_BUILD_PYGMO")
	INSTALL(TARGETS pagmo
		RUNTIME DESTINATION ${LIB_INSTALL_PATH}
		LIBRARY DESTINATION ${LIB_INSTALL_PATH}
	)
ENDIF(BUILD_PYGMO)
