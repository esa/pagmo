ADD_EXECUTABLE(local_torture_test local_torture_test.cpp)
TARGET_LINK_LIBRARIES(local_torture_test ${MANDATORY_LIBRARIES} pagmo_static)
ADD_TEST(local_torture_test local_torture_test)
IF(ENABLE_MPI)
	ADD_EXECUTABLE(mpi_torture_test mpi_torture_test.cpp)
	TARGET_LINK_LIBRARIES(mpi_torture_test ${MANDATORY_LIBRARIES} pagmo_static)
	ADD_TEST(mpi_torture_test ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 10 ${MPIEXEC_PREFLAGS} mpi_torture_test
		${MPIEXEC_POSTFLAGS})
ENDIF(ENABLE_MPI)
