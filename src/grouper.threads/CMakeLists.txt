cmake_minimum_required(VERSION 3.13)

project(grouperThreads${VERSION})

add_subdirectory(testfiles)

#Generate file sources

file(

        GLOB

        source_files

        *

)


#Declare library

add_executable(
        grouperThreads${VERSION}

        ${source_files}

)

target_link_libraries(
        grouperThreads${VERSION}


        SDGlib
        DGElib
        BLRlib
)