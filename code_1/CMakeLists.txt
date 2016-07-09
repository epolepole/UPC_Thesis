#CMakeLists.txt
#set(CMAKE_C_COMPILER /usr/local/berkeley_upc/bin/upcc)

cmake_minimum_required (VERSION 2.6)



#set(UPC "")

if(DEFINED UPC)
	set(CMAKE_C_COMPILER /usr/local/berkeley_upc/bin/upcc)
	message(STATUS "Using Berkeley UPC")
else(DEFINED UPC)
	#set(CMAKE_C_COMPILER /usr/bin/gcc)
	message(STATUS "Using default compiler")
endif(DEFINED UPC)

project (thesis_upc)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build/bin")
set(CMAKE_BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/build")

#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")


if(NOT DEFINED THREADS)
	set(THREADS 4)
	message(STATUS "Setting 4 threads")
endif(NOT DEFINED THREADS)

message(STATUS "Compiling for ${THREADS} threads")

set(CMAKE_C_FLAGS_DEBUG "-T ${THREADS} -g -D__SAVE_ITER__=1")
set(CMAKE_C_FLAGS_RELEASE "-T ${THREADS} -D__SAVE_ITER__=1")


IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -network=mpi")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -network=mpi")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")


include(CheckCSourceCompiles)
include(CheckCSourceRuns)
include(CheckFunctionExists)
include(CheckIncludeFile)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CheckSymbolExists)
include(CheckTypeSize)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

#file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/SetUpData.ini" "${CMAKE_CURRENT_SOURCE_DIR}/Mesh" DESTINATION "${EXECUTABLE_OUTPUT_PATH}")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/Mesh" DESTINATION "${EXECUTABLE_OUTPUT_PATH}")


set(HEAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include/)
include_directories ("./include/")
include_directories(/home/epolepole/Programs/berkeley_upc-2.22.2/upcr_preinclude/ /home/epolepole/Programs/berkeley_upc-2.22.2/)

add_subdirectory (src)

#enable_testing()
#add_test(NAME check_money COMMAND check_money)