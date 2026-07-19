cmake_minimum_required(VERSION 3.18)

#enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

project(DAFS CXX C)
include(CTest)
option(ENABLE_CPLEX "Enable CPLEX backend" OFF)
option(ENABLE_GUROBI "Enable Gurobi backend" OFF)
option(ENABLE_SCIP "Enable SCIP backend" OFF)
option(ENABLE_HIGHS "Enable HiGHS backend" OFF)
option(STATIC_BUILD "Build statically for Linux" OFF)

include_directories("spdlog/include")
include_directories("cxxopts/include")

add_subdirectory(src)
