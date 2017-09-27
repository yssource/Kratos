set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosSwimmingDEMApplication ****")

if(NOT (${DEM_APPLICATION} STREQUAL "ON"))
    message(SEND_ERROR "The DEM Application must be compiled before this one. Set DEM_APPLICATION=ON in the configure file.")
endif(NOT (${DEM_APPLICATION} STREQUAL "ON"))

if(NOT (${FLUID_DYNAMICS_APPLICATION} STREQUAL "ON"))
    message(SEND_ERROR "The Fluid Dynamics Application must be compiled before this one. Set FLUID_DYNAMICS_APPLICATION=ON in the configure file.")
endif(NOT (${FLUID_DYNAMICS_APPLICATION} STREQUAL "ON"))

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

#THIS IS THE ONLY NEW LINE WHEN USING "DERIVED" APPLICATIONS 
#dependecies on other applications
include_directories( ${CMAKE_SOURCE_DIR}/applications/DEM_application )

## generate variables with the sources
set( KRATOS_CUSTOM_DEM_UTILS_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_DEM_utils_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/custom_DEM_utils_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
)


###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosCustomDEMUtilsApplication SHARED ${KRATOS_CUSTOM_DEM_UTILS_APPLICATION_SOURCES})
target_link_libraries(KratosCustomDEMUtilsApplication KratosDEMApplication KratosCore)
set_target_properties(KratosCustomDEMUtilsApplication PROPERTIES PREFIX "")
install(TARGETS KratosCustomDEMUtilsApplication DESTINATION libs )

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosCustomDEMUtilsApplication PROPERTIES SUFFIX .pyd)
    set_target_properties(KratosCustomDEMUtilsApplication PROPERTIES COMPILE_DEFINITIONS "CUSTOM_DEM_UTILS_APPLICATION=EXPORT,API")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosCustomDEMUtilsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosCustomDEMUtilsApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/CustomDEMUtilsApplication.py" DESTINATION KratosMultiphysics )
