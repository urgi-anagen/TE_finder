cmake_minimum_required(VERSION 3.13)

project(SDGlib)
add_subdirectory (testfiles)


#Generate file sources

file(

        GLOB

        source_files

        *.cpp

)


#Declare library

add_library(

        SDGlib

        SHARED

        ${source_files}

)
