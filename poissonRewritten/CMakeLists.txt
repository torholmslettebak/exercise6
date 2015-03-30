project(poisson1d)

cmake_minimum_required(VERSION 2.6)
add_subdirectory(../common common)
# We use the C language
enable_language(C)
SET(CMAKE_COMPILER_IS_GNUCC)
enable_language(Fortran)

add_definitions(${COMMON_DEFINITIONS})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_C_FLAGS}")
include_directories(${COMMON_INCLUDE_PATHS})


INCLUDE(CheckCCompilerFlag)
CHECK_C_COMPILER_FLAG("-std=c99" HAVE_C99)
if (HAVE_C99)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_C_FLAGS} -std=c99 ")
else(HAVE_C99)
  MESSAGE(WARNING "C99 support not available, program may fail to compile!")
ENDIF(HAVE_C99)

include_directories(${COMMON_INCLUDE_PATHS})
# CMake Options to enable/disable MPI and OpenMP
option(ENABLE_OPENMP "Enable OpenMP support?" ON)
option(ENABLE_MPI    "Enable MPI support?" ON)

if(ENABLE_MPI)
  find_package(MPI)
endif()

if(MPI_FOUND)
  add_definitions(-DHAVE_MPI=1)
  set(INCLUDES ${INCLUDES} ${MPI_INCLUDE_PATH})
  set(DEPLIBS_C ${DEPLIBS_C} ${MPI_C_LIBRARIES})
  set(DEPLIBS_F ${DEPLIBS_F} ${MPI_Fortran_LIBRARIES})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_FLAGS}")
endif()

if(ENABLE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    add_definitions(-DHAVE_OPENMP=1)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
  endif()
endif()

# Send the include directories to the compiler
include_directories(${INCLUDES})

add_library(poisson_common poissoncommon.c fst.f)


add_executable(poissonRewrite poissonRewrite.c)
target_link_libraries(poissonRewrite poisson_common ${COMMON_LIBRARIES})