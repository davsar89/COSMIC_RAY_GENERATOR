#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.6 FATAL_ERROR)
project(CR_GEN C CXX Fortran)
enable_language(Fortran)

#execute_process(COMMAND "cd ${PROJECT_SOURCE_DIR}/src/CRY/; make" RESULT_VARIABLE rv)
#message("${PROJECT_SOURCE_DIR}/src/CRY/")
#message("rv='${rv}'")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

#----------------------------------------------------------------------------

find_package(PkgConfig REQUIRED)
find_package(Threads REQUIRED)
pkg_search_module(UUID REQUIRED uuid)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#

include_directories(${PROJECT_SOURCE_DIR}/build/)
include_directories(${PROJECT_SOURCE_DIR}/src/)
include_directories(${PROJECT_SOURCE_DIR}/src/src/)
include_directories(${PROJECT_SOURCE_DIR}/src/src/parma/)
include_directories(${PROJECT_SOURCE_DIR}/src/src/C_functions/)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/src/src/*.hh)
file(GLOB sourcesf90 ${PROJECT_SOURCE_DIR}/src/src/parma2/*.f90)
file(GLOB sourcesC ${PROJECT_SOURCE_DIR}/src/src/C_functions/*.c)
file(GLOB main ${PROJECT_SOURCE_DIR}/src/main.cc)

#${hdf}
#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(CR_GEN ${main} ${sourcesC} ${sourcesf90} ${sources} ${headers})

target_include_directories(CR_GEN PUBLIC ${UUID_INCLUDE_DIRS})
target_link_libraries(CR_GEN ${UUID_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -mfpmath=sse -msse2 -ffp-contract=off")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-512 -fno-range-check")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fp-model=extended")
endif ()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS CR_GEN DESTINATION bin)

