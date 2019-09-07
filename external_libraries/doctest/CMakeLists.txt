set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Generate variables with the sources
set(DOCTEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/doctest.cpp)

###############################################################

add_library(doctest STATIC ${DOCTEST_SOURCES})
target_link_libraries(doctest m)

set(DOCTEST_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(DOCTEST_LIBRARIES doctest)

message(STATUS "doctest : ${DOCTEST_INCLUDE_DIR}" )
