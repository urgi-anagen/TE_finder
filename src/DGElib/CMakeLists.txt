cmake_minimum_required(VERSION 3.5)

project(DGElib)
add_subdirectory (testfiles)


#Generate file sources

file(

        GLOB

        source_files

        *.cpp

)


#Declare library

add_library(

        DGElib

        STATIC

        ${source_files}

)

target_link_libraries(
        DGElib

        SDGlib
)