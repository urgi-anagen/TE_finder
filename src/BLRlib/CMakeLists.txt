cmake_minimum_required(VERSION 3.5)

project(BLRlib)
add_subdirectory (testfiles)


#Generate file sources

file(

        GLOB

        source_files

        *.cpp

)


#Declare library

add_library(
        BLRlib

        SHARED

        ${source_files}
        BLRMatchPath.cpp BLRMatchPath.h BLRMatchJoin.cpp BLRMatchJoin.h)

target_link_libraries(
        BLRlib

        SDGlib
        DGElib
)