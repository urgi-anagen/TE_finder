cmake_minimum_required(VERSION 3.5)

project(blaster${VERSION})


#Generate file sources

file(

        GLOB

        source_files

        *

)


#Declare library

add_executable(
        blaster${VERSION}

        ${source_files}

)

target_link_libraries(
        blaster${VERSION}

        SDGlib
        DGElib
        BLRlib
)