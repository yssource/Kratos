set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

############ module to include sys libs #############
include(InstallRequiredSystemLibraries)

## generate variables with the sources
set( KRATOS_TEST_SOURCES
 	${CMAKE_CURRENT_SOURCE_DIR}/ktests_main.cpp
)

# This needs the new cmake structure to work. Add test libs based on apps if tests building is on
# Comented code below will be used once we implement it:
set( KRATOS_TEST_LIBRARIES
    KratosCoreTests
)

# # Adding kratos application tests from ENV
# if(DEFINED ENV{KRATOS_INT_APP})
#     foreach(APPLICATION_PATH $ENV{KRATOS_INT_APP})
#         add_subdirectory(applications/${APPLICATION_PATH})
#     endforeach(APPLICATION_PATH $ENV{KRATOS_INT_APP})
# endif(DEFINED ENV{KRATOS_INT_APP})

# # Addong kratos application tests from CONFIG
# if(DEFINED KRATOS_INT_APPLICATION)
#     foreach(APPLICATION_PATH ${KRATOS_INT_APPLICATION})
#         add_subdirectory(applications/${APPLICATION_PATH})
#     endforeach(APPLICATION_PATH ${KRATOS_INT_APPLICATION})
# endif(DEFINED KRATOS_INT_APPLICATION)

# # Adding external application tests from ENV
# if(DEFINED ENV{KRATOS_EXT_APP})
#     foreach(APPLICATION_PATH $ENV{KRATOS_EXT_APP})
#         get_filename_component(APPLICATION_NAME ${APPLICATION_PATH} NAME)
#         add_subdirectory(${APPLICATION_PATH} ${CMAKE_CURRENT_BINARY_DIR}/applications/${APPLICATION_NAME})
#     endforeach(APPLICATION_PATH $ENV{KRATOS_EXT_APP})
# endif(DEFINED ENV{KRATOS_EXT_APP})

# # Adding external application tests from CONFIG
# if(DEFINED KRATOS_EXT_APPLICATION)
#     foreach(APPLICATION_PATH ${KRATOS_EXT_APPLICATION})
#         get_filename_component(APPLICATION_NAME ${APPLICATION_PATH} NAME)
#         add_subdirectory(${APPLICATION_PATH} ${CMAKE_CURRENT_BINARY_DIR}/applications/${APPLICATION_NAME})
#     endforeach(APPLICATION_PATH ${KRATOS_EXT_APPLICATION})
# endif(DEFINED KRATOS_EXT_APPLICATION)


###############################################################
add_executable(runtests ${KRATOS_TEST_SOURCES})
target_link_libraries(runtests PUBLIC ${KRATOS_TEST_LIBRARIES})

if(USE_COTIRE MATCHES ON)
    cotire(runtests)
endif(USE_COTIRE MATCHES ON)
install(TARGETS runtests DESTINATION . )