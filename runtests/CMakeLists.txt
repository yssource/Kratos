set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

############ module to include sys libs #############
include(InstallRequiredSystemLibraries)

## generate variables with the sources
set( KRATOS_TEST_SOURCES
 	${CMAKE_CURRENT_SOURCE_DIR}/ktests_main.cpp
)

get_property(KRATOS_TESTING_LIBRARY_LIST GLOBAL PROPERTY KRATOS_TESTING_LIBRARIES)

add_executable(runtests ${KRATOS_TEST_SOURCES})
target_link_libraries(runtests PUBLIC ${KRATOS_TESTING_LIBRARY_LIST})

if(USE_COTIRE MATCHES ON)
    cotire(runtests)
endif(USE_COTIRE MATCHES ON)
install(TARGETS runtests DESTINATION . )