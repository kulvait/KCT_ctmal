cmake_minimum_required(VERSION 3.1)
project(CTIOL)

set(CMAKE_INSTALL_PREFIX $ENV{HOME}/bin)
set(CMAKE_CXX_STANDARD 14)#Supported values are 98, 11 and 14.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
message("Current value of flags is ${CMAKE_CXX_FLAGS} and debug is ${CMAKE_CXX_FLAGS_DEBUG} and release is ${CMAKE_CXX_FLAGS_RELEASE}")
#Default CMAKE_CXX_FLAGS are empty, default for CMAKE_CXX_FLAGS_DEBUG is -g and CMAKE_CXX_FLAGS_RELEASE are not empty
set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_BUILD_TYPE Debug)
set(BUILD_SHARED_LIBS True)

#MKL
find_package(MKL)
include_directories(${MKL_INCLUDE_DIRS})

#Header directories for the project
set(GLOBAL_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
include_directories(${GLOBAL_INCLUDE_DIR})

#Includes relative to submodules directories
include_directories(${CMAKE_SOURCE_DIR}/submodules)

#Plog logger
include_directories(${CMAKE_SOURCE_DIR}/submodules/plog/include)

#Matrix class 
#ssh:  // git@gitlab.stimulate.ovgu.de:2200/robert-frysch/Matrix-Template.git
include_directories(${CMAKE_SOURCE_DIR}/submodules/Matrix-Template)
add_definitions("-DENABLE_FROBENIUS_NORM")#Nor in the frobenius sense

#Projection matrix comparator by Robert Frysch
#To enable Frobenius norm just define ENABLE_FROBENIUS_NORM
#add_definitions("-DENABLE_FROBENIUS_NORM")
#include_directories(${CMAKE_SOURCE_DIR} / submodules / projection - matrix - \
#		     comparator)
#SET(PMC_SRC ${CMAKE_SOURCE_DIR} / submodules / projection - matrix -         \
#     comparator / pmatcachedfilereader.cpp ${CMAKE_SOURCE_DIR} / submodules / \
#	 projection -                                                         \
#     matrix - comparator / pmatcomparator.cpp)
#add_library(pmc ${PMC_SRC})

#CT input output library
FILE( GLOB CTIOL_SRC ${CMAKE_SOURCE_DIR}/submodules/CTIOL/src/*.cpp ${CMAKE_SOURCE_DIR}/submodules/CTIOL/src/ARGPARSE/*.cpp)
add_library(ctiol ${CTIOL_SRC})
include_directories(${CMAKE_SOURCE_DIR}/submodules/CTIOL/include)

#CTMAL library
include_directories(${CMAKE_SOURCE_DIR}/include)
FILE( GLOB CTMAL_SRC ${CMAKE_SOURCE_DIR}/src/*.cpp )
add_library(ctmal ${CTMAL_SRC})
target_link_libraries(ctmal ctiol)
set_target_properties(
    ctmal
    PROPERTIES 
        OUTPUT_NAME "ctmal.so"
        SUFFIX ""
)

add_custom_target(formatWebkit
./formatWebkit
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

#Documentation target
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/doc/doxygen.conf.in)
    set(DOXYGEN_OUT ${CMAKE_SOURCE_DIR}/doc/doxygen.conf)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT})

    # note the option ALL which allows to build the docs together with the application
    add_custom_target( doxygen_doc
        ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating documentation with Doxygen"
        VERBATIM )

    add_custom_target( doxygen
        make
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc/latex
        COMMENT "Generating PDF manual"
        VERBATIM )

        add_dependencies(doxygen doxygen_doc)

else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)

#Testing
include_directories(${CMAKE_SOURCE_DIR}/tests)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CMAKE_SOURCE_DIR}/submodules/Catch2/single_include/catch2)
FILE( GLOB TST_SRC ${CMAKE_SOURCE_DIR}/tests/*.cpp )
add_executable(test_ctmal ${TST_SRC})
target_link_libraries(test_ctmal Catch)
target_link_libraries(test_ctmal ctiol)
target_link_libraries(test_ctmal ctmal)
enable_testing()
add_test(NAME CTMALTest COMMAND test_ctmal)
