cmake_minimum_required(VERSION 3.5)

project(matcher${VERSION})

add_subdirectory(testfiles)

#Generate file sources

file(

        GLOB

        source_files

        *

)


#Declare library

add_executable(
        matcher${VERSION}

        ${source_files}

)

target_link_libraries(
        matcher${VERSION}

        SDGlib
        DGElib
        BLRlib
)