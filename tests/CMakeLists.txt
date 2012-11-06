ADD_EXECUTABLE(local_torture_test local_torture_test.cpp)
TARGET_LINK_LIBRARIES(local_torture_test pagmo_static ${MANDATORY_LIBRARIES})
ADD_TEST(local_torture_test local_torture_test)

ADD_EXECUTABLE(serialization_problems serialization_problems.cpp)
TARGET_LINK_LIBRARIES(serialization_problems pagmo_static ${MANDATORY_LIBRARIES})
ADD_TEST(serialization_problems serialization_problems)

ADD_EXECUTABLE(serialization_algorithms serialization_algorithms.cpp)
TARGET_LINK_LIBRARIES(serialization_algorithms pagmo_static ${MANDATORY_LIBRARIES})
ADD_TEST(serialization_algorithms serialization_algorithms)

IF(ENABLE_MPI)
	ADD_EXECUTABLE(mpi_torture_test mpi_torture_test.cpp)
	TARGET_LINK_LIBRARIES(mpi_torture_test ${MANDATORY_LIBRARIES} pagmo_static)
	ADD_TEST(mpi_torture_test_01 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 10 ${MPIEXEC_PREFLAGS} mpi_torture_test
		${MPIEXEC_POSTFLAGS})
	ADD_TEST(mpi_torture_test_02 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} mpi_torture_test
		${MPIEXEC_POSTFLAGS})
ENDIF(ENABLE_MPI)
