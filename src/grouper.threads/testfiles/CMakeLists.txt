cmake_minimum_required(VERSION 3.5)

project(test_grouperThreads${VERSION})

include_directories(../../DGElib/testfiles
        ../../BLRlib/testfiles)

#Generate file sources

file(

        GLOB

        source_files

        *

)


#Declare library

add_executable(
        test_grouperThreads${VERSION}

        ${source_files}

        ../BLRMemIdxBin.cpp
        ../BLRGroupThreads.cpp
        ../BLRGrouperThreads.cpp
        Test_BLRGroup.cpp Test_BLRGroup.h
        Test_BLRGrouper.cpp Test_BLRGrouper.h)

target_link_libraries(
        test_grouperThreads${VERSION}

        ${CPPUNIT_LIBRARY}
        SDGlib
        DGElib
        BLRlib
)