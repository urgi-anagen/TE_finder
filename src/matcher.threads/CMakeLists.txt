cmake_minimum_required(VERSION 3.5)

project(matcherThreads${VERSION})

add_subdirectory(testfiles)

#Generate file sources

file(

        GLOB

        source_files

        *

)


#Declare library

add_executable(
        matcherThreads${VERSION}

        ${source_files}

)

target_link_libraries(
        matcherThreads${VERSION}


        SDGlib
        DGElib
        BLRlib
)

install (TARGETS matcherThreads${VERSION} DESTINATION ${CMAKE_SOURCE_DIR}/bin)