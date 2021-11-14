cmake_minimum_required(VERSION 3.8)

#enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

project(DAFS CXX C)

include_directories("spdlog/include")
include_directories("cxxopts/include")

add_subdirectory(src)