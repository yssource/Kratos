set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMGISApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/MGISApplication )


## Define application core sources (pure C++ parts)
set( KRATOS_MGIS_APPLICATION_CORE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/mgis_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mgis_application_variables.cpp

  # CUSTOM UTILITIES

  # CUSTOM PROCESSES

)

## generate variables with the testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
  file(GLOB_RECURSE KRATOS_MGIS_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## Define python interface sources
set( KRATOS_MGIS_APPLICATION_PYTHON_INTERFACE_SOURCES
  ## CUSTOM PYTHON
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/mgis_python_application.cpp
)

###############################################################################
## MGISApplication core library (C++ parts)

# If build mode is not defined, assume Release
if(NOT MGIS_INCLUDE_DIR)
  message(FATAL_ERROR "Undefined MGIS_INCLUDE_DIR")
endif(NOT MGIS_INCLUDE_DIR)
if(NOT MGIS_LIBRARY)
  message(FATAL_ERROR "Undefined MGIS_LIBRARY")
endif(NOT MGIS_LIBRARY)
add_library(KratosMGISCore SHARED ${KRATOS_MGIS_APPLICATION_CORE_SOURCES} ${KRATOS_MGIS_APPLICATION_TESTING_SOURCES})
INCLUDE_DIRECTORIES(${MGIS_INCLUDE_DIR})
target_link_libraries(KratosMGISCore PUBLIC KratosCore KratosStructuralMechanicsCore ${MGIS_LIBRARY})
set_target_properties(KratosMGISCore PROPERTIES COMPILE_DEFINITIONS "MGIS_APPLICATION=EXPORT,API")

## MGISApplication python module
pybind11_add_module(KratosMGISApplication MODULE THIN_LTO ${KRATOS_MGIS_APPLICATION_PYTHON_INTERFACE_SOURCES})
target_link_libraries(KratosMGISApplication PUBLIC KratosMGISCore)
set_target_properties(KratosMGISApplication PROPERTIES PREFIX "")

###############################################################################
# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosMGISApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosMGISApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN "*.csv")
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
    PATTERN "*.git" EXCLUDE
    PATTERN "*.c" EXCLUDE
    PATTERN "*.h" EXCLUDE
    PATTERN "*.cpp" EXCLUDE
    PATTERN "*.hpp" EXCLUDE
  )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

###############################################################################
if(USE_COTIRE MATCHES ON)
    cotire(KratosMGISCore)
    cotire(KratosMGISApplication)
endif(USE_COTIRE MATCHES ON)

###############################################################################
## installing the resulting libraries

install(TARGETS KratosMGISCore DESTINATION libs )
install(TARGETS KratosMGISApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MGISApplication.py" DESTINATION "KratosMultiphysics/MGISApplication" RENAME "__init__.py")
