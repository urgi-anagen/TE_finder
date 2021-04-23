cmake_minimum_required(VERSION 3.5)

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

        STATIC

        ${source_files}

        BioSeq.cpp BioSeq.h
        FastaIstream.cpp FastaIstream.h
        FastaOstream.cpp FastaOstream.h)
