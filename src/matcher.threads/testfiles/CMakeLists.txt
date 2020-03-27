cmake_minimum_required(VERSION 3.5)

project(test_matcherThreads${VERSION})

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
        test_matcherThreads${VERSION}

        ${source_files}

        ../BLRMatcherThreads.cpp
        ../BLRMatcherThreads.h
        ../BLRMatcherThreadsParameter.cpp
        ../BLRMatcherThreadsParameter.h
        Test_matcherThreads.cpp Test_matcherThreads.h
)

target_link_libraries(
        test_matcherThreads${VERSION}

        ${CPPUNIT_LIBRARY}
        SDGlib
        DGElib
        BLRlib
)