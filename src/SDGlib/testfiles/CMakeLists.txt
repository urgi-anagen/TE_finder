cmake_minimum_required(VERSION 3.5)

project(SDGlib_test)


#Generate file sources

file(

        GLOB

        source_files

        *.cpp

)


#Declare library

add_executable(
        test_SDG

        ${source_files}

)

target_link_libraries(
        test_SDG

        ${CPPUNIT_LIBRARY}
        SDGlib
)