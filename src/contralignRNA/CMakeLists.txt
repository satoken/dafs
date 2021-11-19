add_library(contralignRNA
    STATIC
    ../contralign/wrapper.cpp ../contralign/MultiSequence.cpp ../contralign/Sequence.cpp ../contralign/Utilities.cpp 
)
add_compile_definitions(RNA=1)