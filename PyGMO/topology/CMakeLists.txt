ADD_LIBRARY(_topology SHARED topology.cpp)
PYTHON_CPPMODULE(_topology)
INSTALL(FILES __init__.py DESTINATION ${PYGMO_INSTALL_PATH}/topology)
