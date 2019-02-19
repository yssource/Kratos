set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosRANSConstitutiveLawsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_RANS_CONSTITUTIVE_LAWS_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/rans_constitutive_laws_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/rans_constitutive_laws_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/rans_constitutive_laws_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/rans_constitutive_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/evm_k_epsilon/evm_k_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/evm_k_epsilon/evm_epsilon_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/evm_k_epsilon_wall_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/turbulence_eddy_viscosity_model_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
)

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosRANSConstitutiveLawsApplication MODULE THIN_LTO ${KRATOS_RANS_CONSTITUTIVE_LAWS_APPLICATION_SOURCES})
target_link_libraries(KratosRANSConstitutiveLawsApplication PRIVATE KratosCore )
set_target_properties(KratosRANSConstitutiveLawsApplication PROPERTIES PREFIX "")

if(USE_COTIRE MATCHES ON)
    cotire(KratosRANSConstitutiveLawsApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosRANSConstitutiveLawsApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosRANSConstitutiveLawsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosRANSConstitutiveLawsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/RANSConstitutiveLawsApplication.py" DESTINATION "KratosMultiphysics/RANSConstitutiveLawsApplication" RENAME "__init__.py")
