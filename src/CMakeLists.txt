# PaGMO lib source files.
SET(PAGMO_LIB_SRC_LIST
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/base.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/ihs.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/monte_carlo.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/algorithm/null.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/migration/base.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/migration/base_r_policy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/migration/base_s_policy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/migration/best_s_policy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/migration/fair_r_policy.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/base.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/golomb_ruler.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/himmelblau.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/knapsack.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/paraboloid.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/problem/rastrigin.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/barabasi_albert.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/base.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/custom.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/erdos_renyi.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/fully_connected.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/one_way_ring.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/ring.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/unconnected.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topology/watts_strogatz.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/archipelago.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/island.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/population.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/rng.cpp
)

# Boost  thread library.
IF(UNIX)
	SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST}
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/pthread/once.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/pthread/thread.cpp
	)
ELSE(UNIX)
	SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST}
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/win32/thread.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/win32/tss_dll.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/win32/tss_pe.cpp
	)
ENDIF(UNIX)

# Boost serialization library.
SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST}
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/archive_exception.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/basic_serializer_map.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/binary_oarchive.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/polymorphic_oarchive.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/void_cast.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/basic_archive.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/basic_text_iprimitive.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/binary_wiarchive.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/shared_ptr_helper.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/xml_archive_exception.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/basic_iarchive.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/basic_text_oprimitive.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/binary_woarchive.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/stl_port.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/xml_grammar.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/basic_iserializer.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/text_iarchive.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/xml_iarchive.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/basic_oarchive.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/extended_type_info.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/text_oarchive.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/xml_oarchive.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/basic_oserializer.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/basic_xml_archive.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/extended_type_info_no_rtti.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/basic_pointer_iserializer.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/basic_xml_grammar.ipp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/extended_type_info_typeid.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/basic_pointer_oserializer.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/binary_iarchive.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/polymorphic_iarchive.cpp
)

# Disable UTF-8 serialization in MinGW, since the latest stable version does not support it.
IF(NOT MINGW)
	SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST}
		${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/basic_text_wiprimitive.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/basic_text_woprimitive.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/text_wiarchive.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/text_woarchive.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/utf8_codecvt_facet.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/xml_wgrammar.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/xml_wiarchive.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/xml_woarchive.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/codecvt_null.cpp
	)
ENDIF(NOT MINGW)

# Create a pagmo_static library if main is requested.
IF(BUILD_MAIN)
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
	SET(ADDITIONAL_FLAGS "-DBOOST_PYTHON_SOURCE -DBOOST_THREAD_BUILD_DLL -DBOOST_SERIALIZATION_DYN_LINK=1 -DBOOST_WSERIALIZATION_DYN_LINK=1 -DPAGMO_BUILD_PYGMO -DPAGMO_DLL_EXPORT_API")
	SET_TARGET_PROPERTIES(pagmo PROPERTIES COMPILE_FLAGS "${ADDITIONAL_FLAGS}")
	
	INSTALL(TARGETS pagmo
		RUNTIME DESTINATION ${LIB_INSTALL_PATH}
		LIBRARY DESTINATION ${LIB_INSTALL_PATH}
	)
ENDIF(BUILD_PYGMO)
