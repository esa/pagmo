ADD_EXECUTABLE(best_solutions_test best_solutions_test.cpp)
TARGET_LINK_LIBRARIES(best_solutions_test ${MANDATORY_LIBRARIES} pagmo_static)
ADD_TEST(best_solutions_test best_solutions_test)

ADD_EXECUTABLE(local_torture_test local_torture_test.cpp)
TARGET_LINK_LIBRARIES(local_torture_test ${MANDATORY_LIBRARIES} pagmo_static)
ADD_TEST(local_torture_test local_torture_test)

ADD_EXECUTABLE(serialization_problems serialization_problems.cpp)
TARGET_LINK_LIBRARIES(serialization_problems ${MANDATORY_LIBRARIES} pagmo_static)
ADD_TEST(serialization_problems serialization_problems)

ADD_EXECUTABLE(serialization_algorithms serialization_algorithms.cpp)
TARGET_LINK_LIBRARIES(serialization_algorithms ${MANDATORY_LIBRARIES} pagmo_static)
ADD_TEST(serialization_algorithms serialization_algorithms)

ADD_EXECUTABLE(test_shifted test_shifted.cpp)
TARGET_LINK_LIBRARIES(test_shifted ${MANDATORY_LIBRARIES} pagmo_static)
ADD_TEST(test_shifted test_shifted)

ADD_EXECUTABLE(test_rotated test_rotated.cpp)
TARGET_LINK_LIBRARIES(test_rotated ${MANDATORY_LIBRARIES} pagmo_static)
ADD_TEST(test_rotated test_rotated)

ADD_EXECUTABLE(test_noisy test_noisy.cpp)
TARGET_LINK_LIBRARIES(test_noisy ${MANDATORY_LIBRARIES} pagmo_static)
ADD_TEST(test_noisy test_noisy)

ADD_EXECUTABLE(hypervolume_test hypervolume_test.cpp)
TARGET_LINK_LIBRARIES(hypervolume_test ${MANDATORY_LIBRARIES} pagmo_static)
ADD_TEST(hypervolume_test hypervolume_test)

ADD_EXECUTABLE(serialization_hypervolume serialization_hypervolume.cpp)
TARGET_LINK_LIBRARIES(serialization_hypervolume ${MANDATORY_LIBRARIES} pagmo_static)
ADD_TEST(serialization_hypervolume serialization_hypervolume)

ADD_EXECUTABLE(test_robust test_robust.cpp)
TARGET_LINK_LIBRARIES(test_robust ${MANDATORY_LIBRARIES} pagmo_static)
ADD_TEST(test_robust test_robust)

ADD_EXECUTABLE(test_racing test_racing.cpp)
TARGET_LINK_LIBRARIES(test_racing ${MANDATORY_LIBRARIES} pagmo_static)
ADD_TEST(test_racing test_racing)

ADD_EXECUTABLE(test_racing_algorithm test_racing_algorithm.cpp)
TARGET_LINK_LIBRARIES(test_racing_algorithm ${MANDATORY_LIBRARIES} pagmo_static)
ADD_TEST(test_racing_algorithm test_racing_algorithm)

ADD_EXECUTABLE(test_tsp test_tsp.cpp)
TARGET_LINK_LIBRARIES(test_tsp ${MANDATORY_LIBRARIES} pagmo_static)
ADD_TEST(test_tsp test_tsp)

ADD_EXECUTABLE(test_archipelago test_archipelago.cpp)
TARGET_LINK_LIBRARIES(test_archipelago ${MANDATORY_LIBRARIES} pagmo_static)
ADD_TEST(test_archipelago test_archipelago)

IF(ENABLE_MPI)
	ADD_EXECUTABLE(mpi_torture_test mpi_torture_test.cpp)
        TARGET_LINK_LIBRARIES(mpi_torture_test ${MANDATORY_LIBRARIES} pagmo_static)
	ADD_TEST(mpi_torture_test_01 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 10 ${MPIEXEC_PREFLAGS} ./mpi_torture_test
		${MPIEXEC_POSTFLAGS})
	ADD_TEST(mpi_torture_test_02 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} ./mpi_torture_test
		${MPIEXEC_POSTFLAGS})
ENDIF(ENABLE_MPI)
