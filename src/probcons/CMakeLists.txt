add_library(probcons
    STATIC
    wrapper.cpp
)
add_compile_definitions(NumInsertStates=2 VERSION="1.12")