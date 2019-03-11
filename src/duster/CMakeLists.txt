cmake_minimum_required(VERSION 3.13)

project(duster${VERSION})

add_subdirectory(testfiles)

#Generate file sources

file(

        GLOB

        source_files

        *

)


#Declare library

add_executable(
        duster${VERSION}

        ${source_files}

)

target_link_libraries(
        duster${VERSION}

        SDGlib
        DGElib
        BLRlib
)