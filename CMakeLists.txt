cmake_minimum_required(VERSION 3.5)
project(TE_finder)

set(CMAKE_CXX_STANDARD 11)
set (VERSION "2.29")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERSION=${VERSION} -pthread -Wall")

# Search for cppunit library

FIND_PATH(CPPUNIT_INCLUDE_DIR cppunit/Test.h)
FIND_LIBRARY(CPPUNIT_LIBRARY NAMES cppunit)

IF (CPPUNIT_INCLUDE_DIR AND CPPUNIT_LIBRARY)
    SET(CPPUNIT_FOUND TRUE)
ENDIF (CPPUNIT_INCLUDE_DIR AND CPPUNIT_LIBRARY)

IF (CPPUNIT_FOUND)

    # show which CppUnit was found only if not quiet
    IF (NOT CppUnit_FIND_QUIETLY)
        MESSAGE(STATUS "Found CppUnit: ${CPPUNIT_LIBRARY}")
    ENDIF (NOT CppUnit_FIND_QUIETLY)

ELSE (CPPUNIT_FOUND)

    # fatal error if CppUnit is required but not found
    IF (CppUnit_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find CppUnit")
    ENDIF (CppUnit_FIND_REQUIRED)

ENDIF (CPPUNIT_FOUND)

#-----

include_directories(
        ${CPPUNIT_INCLUDE_DIR}
        src/SDGlib
        src/DGElib
        src/BLRlib)

link_directories(
        ${CPPUNIT_LIBRARY_DIR}
        src/SDGlib
        src/DGElib
        src/BLRlib)


add_subdirectory (src/SDGlib)
add_subdirectory (src/DGElib)
add_subdirectory (src/BLRlib)
add_subdirectory (src/blaster)
add_subdirectory (src/matcher)
add_subdirectory (src/grouper.threads)
add_subdirectory (src/duster)
add_subdirectory (src/tools)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)