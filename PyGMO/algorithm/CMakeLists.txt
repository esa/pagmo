ADD_LIBRARY(_algorithm SHARED algorithm.cpp)
PYTHON_CPPMODULE(_algorithm)
INSTALL(FILES __init__.py DESTINATION ${PYGMO_INSTALL_PATH}/algorithm)
