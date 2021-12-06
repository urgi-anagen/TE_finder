cmake_minimum_required(VERSION 3.5)

project(test_matcher${VERSION})


#Generate file sources

file(

        GLOB

        source_files

        *

)


#Declare library

add_executable(
        test_matcher${VERSION}

        ${source_files}

)

target_link_libraries(
        test_matcher${VERSION}

        ${CPPUNIT_LIBRARY}
        SDGlib
        DGElib
        BLRlib
)