cmake_minimum_required(VERSION 3.5)

project(test_Tools)


#Generate file sources

file(

        GLOB

        source_files

        *.cpp

)


#Declare library

add_executable(
        test_Tools

        ${source_files}
        test_Tools.cpp
        test_Tools.h
)

target_link_libraries(
        test_Tools

        ${CPPUNIT_LIBRARY}
        SDGlib
        DGElib
        BLRlib
)