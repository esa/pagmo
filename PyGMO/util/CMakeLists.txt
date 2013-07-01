ADD_LIBRARY(_hypervolume SHARED hypervolume.cpp)
PYTHON_CPPMODULE(_hypervolume)
INSTALL(FILES __init__.py DESTINATION ${PYGMO_INSTALL_PATH}/hypervolume)
