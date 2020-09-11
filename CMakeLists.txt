cmake_minimum_required(VERSION 3.8)

#enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

project(DAFS CXX C)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
add_definitions(-DPACKAGE="dafs")
add_definitions(-DPACKAGE_VERSION="0.0.4")
add_definitions(-DVERSION="0.0.4")
add_definitions(-DRNA=1)
add_definitions(-DNumInsertStates=1)

find_package(PkgConfig REQUIRED)

pkg_check_modules(VIENNARNA REQUIRED RNAlib2)
include_directories(${VIENNARNA_INCLUDE_DIRS})
add_definitions(-DHAVE_VIENNA18)
add_definitions(-DHAVE_VIENNA20)

if(ENABLE_CPLEX)
  find_package(CPLEX)
  if (CPLEX_FOUND)
    add_definitions(-DWITH_CPLEX)
    include_directories(${CPLEX_INCLUDE_DIR})
  endif()
elseif(ENABLE_GUROBI)
  find_package(GUROBI)
  if (GUROBI_FOUND)
    add_definitions(-DWITH_GUROBI)
    include_directories(${GUROBI_INCLUDE_DIR})
  endif()
else()
  find_package(GLPK 4.33)
  if(GLPK_FOUND)
    find_library(ZLIB_LIBRARY z)
    find_library(GMP_LIBRARY gmp)
    find_library(DL_LIBRARY dl)
    find_library(LTDL_LIBRARY ltdl)
    find_library(MPFR_LIBRARY mpfr)
    add_definitions(-DWITH_GLPK)
    include_directories(${GLPK_INCLUDE_DIR})
  endif()
endif()

add_executable(
    dafs
    src/alifold.cpp src/align.cpp src/boltzmann_param.c src/cmdline.c src/dafs.cpp src/fa.cpp src/fold.cpp
    src/ip.cpp src/ipknot.cpp src/needleman_wunsch.cpp src/nussinov.cpp src/tree.cpp
    src/contrafold/SStruct.cpp src/contrafold/Utilities.cpp src/contrafold/wrapper.cpp
    src/contralign/MultiSequence.cpp src/contralign/Sequence.cpp src/contralign/Utilities.cpp src/contralign/wrapper.cpp
    src/probconsRNA/wrapper.cpp
)
target_link_libraries(dafs 
  ${VIENNARNA_LDFLAGS} ${GLPK_LIBRARY} ${ZLIB_LIBRARY} ${GMP_LIBRARY} ${DL_LIBRARY} ${LTDL_LIBRARY} ${MPFR_LIBRARY}
  ${CPLEX_LIBRARY} ${GUROBI_LIBRARY} ${GUROBI_CPP_LIBRARY})
target_compile_features(dafs PUBLIC cxx_std_17)

install(TARGETS dafs RUNTIME DESTINATION bin)