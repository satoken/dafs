add_subdirectory(contrafold)
add_subdirectory(probconsRNA)
add_subdirectory(contralign)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
add_definitions(-DPACKAGE="dafs")
add_definitions(-DPACKAGE_VERSION="0.0.4")
add_definitions(-DVERSION="0.0.4")

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
elseif(ENABLE_GLPK)
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

find_package(Boost 1.40 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(
    dafs
    alifold.cpp align.cpp boltzmann_param.c cmdline.c dafs.cpp fa.cpp fold.cpp
    ip.cpp ipknot.cpp needleman_wunsch.cpp nussinov.cpp tree.cpp
    alignment_graph.cpp gradient_method.cpp graph_library.cpp
)
target_link_libraries(dafs contrafold probcons contralign
  ${VIENNARNA_LDFLAGS} ${GLPK_LIBRARY} ${ZLIB_LIBRARY} ${GMP_LIBRARY} ${DL_LIBRARY} ${LTDL_LIBRARY} ${MPFR_LIBRARY}
  ${CPLEX_LIBRARY} ${GUROBI_LIBRARY} ${GUROBI_CPP_LIBRARY} ${CMAKE_DL_LIBS})
target_compile_features(dafs PUBLIC cxx_std_17)

install(TARGETS dafs RUNTIME DESTINATION bin)
