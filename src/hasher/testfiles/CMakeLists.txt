cmake_minimum_required(VERSION 3.5)

project(test_hasher${VERSION})

include_directories(..)
#Generate file sources

file(

        GLOB

        source_files

        *
        ../Hasher.cpp
        ../Hasher.h
)


#Declare library

add_executable(
        test_hasher${VERSION}

        ${source_files}
)

target_link_libraries(
        test_hasher${VERSION}

        ${CPPUNIT_LIBRARY}
        SDGlib
        DGElib
        BLRlib
)