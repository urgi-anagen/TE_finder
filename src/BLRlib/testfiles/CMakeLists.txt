cmake_minimum_required(VERSION 3.5)

project(BLRlib_test)


#Inclusion des en-têtes publics

include_directories(.
        ..
        ../../SDGlib
        ../../DGElib)

#Generate file sources

file(

        GLOB

        source_files

        *.cpp

)


#Declare library

add_executable(
        test_BLR

        ${source_files}

        Test_BLRUtils.cpp Test_BLRUtils.h)

target_link_libraries(
        test_BLR

        ${CPPUNIT_LIBRARY}
        SDGlib
        DGElib
        BLRlib
)