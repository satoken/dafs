add_subdirectory(probconsRNA)
add_subdirectory(probcons)
add_subdirectory(contralignRNA)
add_subdirectory(contralign)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
add_compile_definitions(PACKAGE="dafs" PACKAGE_VERSION="0.0.4" VERSION="0.0.4")

if(ENABLE_CPLEX)
  find_package(CPLEX)
  if (CPLEX_FOUND)
    add_definitions(-DWITH_CPLEX)
    include_directories(${CPLEX_INCLUDE_DIR})
    find_package(Threads REQUIRED)
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
    align.cpp dafs.cpp fa.cpp
    ip.cpp needleman_wunsch.cpp tree.cpp
    alignment_graph.cpp gradient_method.cpp graph_library.cpp
)
target_link_libraries(dafs probconsRNA probcons contralign contralignRNA
  ${GLPK_LIBRARY} ${ZLIB_LIBRARY} ${GMP_LIBRARY} ${DL_LIBRARY} ${LTDL_LIBRARY} ${MPFR_LIBRARY}
  ${CPLEX_LIBRARY} ${GUROBI_LIBRARY} ${GUROBI_CPP_LIBRARY} ${CMAKE_THREAD_LIBS_INIT} ${CMAKE_DL_LIBS})
target_compile_features(dafs PUBLIC cxx_std_17)

install(TARGETS dafs RUNTIME DESTINATION bin)