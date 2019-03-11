cmake_minimum_required(VERSION 3.13)

project(TE_finder)


set(CMAKE_CXX_STANDARD 11)
set (VERSION "2.28")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERSION=${VERSION} -pthread -Wall")

#set cppunit dirs
#set(CPPUNIT_INCLUDE_DIR "/usr/local/Cellar/cppunit/1.14.0/include")
#set(CPPUNIT_LIBRARY DIR "/usr/local/Cellar/cppunit/1.14.0/lib")

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
        SDGlib
        DGElib
        BLRlib)

link_directories(
        ${CPPUNIT_LIBRARY_DIR}
        SDGlib
        DGElib
        BLRlib)


add_subdirectory (SDGlib)
add_subdirectory (DGElib)
add_subdirectory (BLRlib)
add_subdirectory (blaster)
add_subdirectory (matcher)
add_subdirectory (matcher.threads)
add_subdirectory (grouper.threads)
add_subdirectory (duster)
add_subdirectory (hasher2)
add_subdirectory (tools)

# add the install targets
set(CMAKE_INSTALL_PREFIX ../ CACHE PATH "Install here" FORCE)

install (TARGETS blaster${VERSION} DESTINATION bin)
install (TARGETS matcher.threads${VERSION} DESTINATION bin)
install (TARGETS grouperThreads${VERSION} DESTINATION bin)
install (TARGETS duster${VERSION} DESTINATION ../bin)

#include(CTest)

# does test
#add_test (SDG test_SDG "Ok")

# build a CPack driven installer package
# run: cpack --config CPackConfig.cmake
#include (InstallRequiredSystemLibraries)
#set (CPACK_RESOURCE_FILE_LICENSE
#        "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
#set (CPACK_PACKAGE_VERSION_MAJOR "${Tutorial_VERSION_MAJOR}")
#set (CPACK_PACKAGE_VERSION_MINOR "${Tutorial_VERSION_MINOR}")
#include (CPack)