cmake_minimum_required(VERSION 3.13)

project(test_duster${VERSION})

include_directories(..)
#Generate file sources

file(

        GLOB

        source_files

        *

)


#Declare library

add_executable(
        test_duster${VERSION}

        ${source_files}
        ../Duster.cpp
)

target_link_libraries(
        test_duster${VERSION}

        ${CPPUNIT_LIBRARY}
        SDGlib
        DGElib
        BLRlib
)