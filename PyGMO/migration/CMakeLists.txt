ADD_LIBRARY(_migration SHARED migration.cpp)
PYTHON_CPPMODULE(_migration)
INSTALL(FILES __init__.py DESTINATION ${PYGMO_INSTALL_PATH}/migration)
