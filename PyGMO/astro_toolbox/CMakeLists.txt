ADD_LIBRARY(_astro_toolbox SHARED astro_toolbox.cpp)
PYTHON_CPPMODULE(_astro_toolbox)
INSTALL(FILES __init__.py DESTINATION ${PYGMO_INSTALL_PATH}/astro_toolbox)
