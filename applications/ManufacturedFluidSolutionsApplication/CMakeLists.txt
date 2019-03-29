set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosManufacturedFluidSolutionsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_MANUFACTURED_FLUID_SOLUTIONS_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/manufactured_fluid_solutions_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/manufactured_fluid_solutions_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/manufactured_fluid_solutions_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_manufactured_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_manufactured/manufactured_solution.cpp
)

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosManufacturedFluidSolutionsApplication MODULE THIN_LTO ${KRATOS_MANUFACTURED_FLUID_SOLUTIONS_APPLICATION_SOURCES})
target_link_libraries(KratosManufacturedFluidSolutionsApplication PRIVATE KratosCore )
set_target_properties(KratosManufacturedFluidSolutionsApplication PROPERTIES PREFIX "")

if(USE_COTIRE MATCHES ON)
    cotire(KratosManufacturedFluidSolutionsApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosManufacturedFluidSolutionsApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosManufacturedFluidSolutionsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosManufacturedFluidSolutionsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ManufacturedFluidSolutionsApplication.py" DESTINATION "KratosMultiphysics/ManufacturedFluidSolutionsApplication" RENAME "__init__.py")
