ADD_EXECUTABLE(hm_2_asteroids hm_2_asteroids.cpp)
TARGET_LINK_LIBRARIES(hm_2_asteroids ${MANDATORY_LIBRARIES} pagmo_static)

ADD_EXECUTABLE(migrate_or_not migrate_or_not.cpp)
TARGET_LINK_LIBRARIES(migrate_or_not ${MANDATORY_LIBRARIES} pagmo_static)

IF(ENABLE_SNOPT)
	ADD_EXECUTABLE(gtoc_2_turin gtoc_2_turin.cpp)
	TARGET_LINK_LIBRARIES(gtoc_2_turin ${MANDATORY_LIBRARIES} pagmo_static)
ENDIF(ENABLE_SNOPT)


