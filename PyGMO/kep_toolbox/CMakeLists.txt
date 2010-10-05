ADD_LIBRARY(_kep_toolbox SHARED kep_toolbox.cpp)
PYTHON_CPPMODULE(_kep_toolbox)
INSTALL(FILES __init__.py DESTINATION ${PYGMO_INSTALL_PATH}/kep_toolbox)
