add_library(contralign
    STATIC
    wrapper.cpp MultiSequence.cpp Sequence.cpp Utilities.cpp 
)
add_compile_definitions(RNA=1)