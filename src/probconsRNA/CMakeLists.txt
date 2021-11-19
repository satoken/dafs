add_library(probconsRNA
    STATIC
    wrapper.cpp
)
add_compile_definitions(NumInsertStates=1 VERSION="1.10")