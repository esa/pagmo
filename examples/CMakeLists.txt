IF(ENABLE_GTOP_DATABASE)

ADD_EXECUTABLE(hm_2_asteroids hm_2_asteroids.cpp)
TARGET_LINK_LIBRARIES(hm_2_asteroids pagmo_static ${MANDATORY_LIBRARIES})

ADD_EXECUTABLE(migrate_or_not migrate_or_not.cpp)
TARGET_LINK_LIBRARIES(migrate_or_not pagmo_static ${MANDATORY_LIBRARIES})

IF(ENABLE_SNOPT)
	ADD_EXECUTABLE(gtoc_2_turin gtoc_2_turin.cpp)
        TARGET_LINK_LIBRARIES(gtoc_2_turin pagmo_static ${MANDATORY_LIBRARIES})
ENDIF(ENABLE_SNOPT)

ENDIF(ENABLE_GTOP_DATABASE)

IF(ENABLE_GSL)
	ADD_EXECUTABLE(evolve_spheres evolve_spheres.cpp)
        TARGET_LINK_LIBRARIES(evolve_spheres pagmo_static ${MANDATORY_LIBRARIES})
	ADD_EXECUTABLE(evolve_spheres_racing evolve_spheres_racing.cpp)
        TARGET_LINK_LIBRARIES(evolve_spheres_racing pagmo_static ${MANDATORY_LIBRARIES})
ENDIF(ENABLE_GSL)

ADD_EXECUTABLE(constraints_handling constraints_handling.cpp)
TARGET_LINK_LIBRARIES(constraints_handling pagmo_static ${MANDATORY_LIBRARIES})

ADD_EXECUTABLE(cstrs_death_penalty cstrs_death_penalty.cpp)
TARGET_LINK_LIBRARIES(cstrs_death_penalty pagmo_static ${MANDATORY_LIBRARIES})

ADD_EXECUTABLE(cstrs_co_evolution cstrs_co_evolution.cpp)
TARGET_LINK_LIBRARIES(cstrs_co_evolution pagmo_static ${MANDATORY_LIBRARIES})

ADD_EXECUTABLE(cstrs_self_adaptive cstrs_self_adaptive.cpp)
TARGET_LINK_LIBRARIES(cstrs_self_adaptive pagmo_static ${MANDATORY_LIBRARIES})

ADD_EXECUTABLE(cstrs_self_adaptive_island cstrs_self_adaptive_island.cpp)
TARGET_LINK_LIBRARIES(cstrs_self_adaptive_island pagmo_static ${MANDATORY_LIBRARIES})

ADD_EXECUTABLE(cstrs_immune_system cstrs_immune_system.cpp)
TARGET_LINK_LIBRARIES(cstrs_immune_system pagmo_static ${MANDATORY_LIBRARIES})

ADD_EXECUTABLE(cstrs_core cstrs_core.cpp)
TARGET_LINK_LIBRARIES(cstrs_core pagmo_static ${MANDATORY_LIBRARIES})
