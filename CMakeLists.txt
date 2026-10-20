cmake_minimum_required(VERSION 3.18)

#enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

project(DAFS CXX C)
include(CTest)
option(ENABLE_ILP "Discover and link an integer-programming solver" ON)
option(ENABLE_CPLEX "Enable CPLEX backend" OFF)
option(ENABLE_GUROBI "Enable Gurobi backend" OFF)
option(ENABLE_SCIP "Enable SCIP backend" OFF)
option(ENABLE_HIGHS "Enable HiGHS backend" OFF)
option(STATIC_BUILD "Build statically for Linux" OFF)

include_directories("spdlog/include")
include_directories("cxxopts/include")

# The bundled snapshot needs no external IPknot checkout or probability engines.
# An empty source path selects an installed IPknotDecoder package instead.
set(IPKNOT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/vendor/ipknot" CACHE PATH
    "IPknot checkout containing decoder/ (empty: find installed IPknotDecoder)")
if(IPKNOT_SOURCE_DIR)
  if(NOT EXISTS "${IPKNOT_SOURCE_DIR}/decoder/CMakeLists.txt")
    message(FATAL_ERROR
      "IPknot decoder not found in ${IPKNOT_SOURCE_DIR}. Set IPKNOT_SOURCE_DIR "
      "to an IPknot dev checkout, or leave it empty and set CMAKE_PREFIX_PATH "
      "to an installed IPknotDecoder package.")
  endif()
  add_subdirectory("${IPKNOT_SOURCE_DIR}/decoder" ipknot-decoder)
else()
  find_package(IPknotDecoder 1.1 CONFIG REQUIRED)
endif()

add_subdirectory(src)
