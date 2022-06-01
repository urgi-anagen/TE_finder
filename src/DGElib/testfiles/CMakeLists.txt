cmake_minimum_required(VERSION 3.5)

project(DGElib_test)


#Inclusion des en-têtes publics


#Generate file sources

file(

        GLOB

        source_files

        *.cpp

)

#Declare library

add_executable(
        test_DGE

        ${source_files}

        Test_Graph.cpp Test_Graph.h Test_FastExtAlign.cpp Test_FastExtAlign.h)

target_link_libraries(
        test_DGE

        ${CPPUNIT_LIBRARY}
        SDGlib
        DGElib
)