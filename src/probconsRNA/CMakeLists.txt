add_library(probcons
    STATIC
    wrapper.cpp
)
add_definitions(-DNumInsertStates=1)
add_definitions(-DVERSION="1.10")