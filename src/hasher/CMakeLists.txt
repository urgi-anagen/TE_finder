cmake_minimum_required(VERSION 3.5)

project(hasher${VERSION})

add_subdirectory(testfiles)

#Generate file sources

file(

        GLOB

        source_files

        *

)


#Declare library

add_executable(
        hasher${VERSION}

        ${source_files}

        ../duster/Duster.cpp
        ../duster/Duster.h
)

target_link_libraries(
        hasher${VERSION}

        SDGlib
        DGElib
        BLRlib
)

install (TARGETS hasher${VERSION} DESTINATION ${CMAKE_SOURCE_DIR}/bin)