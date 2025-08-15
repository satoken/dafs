# CMakeLists.txt for LinearAlign library

# Set source files for LinearAlign with simplified BeamAlign implementation
set(LINEARALIGN_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/BeamAlignSimple.cpp
)

# Create static library
add_library(linearalign STATIC ${LINEARALIGN_SOURCES})

# Include directories
target_include_directories(linearalign PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Set C++ standard
target_compile_features(linearalign PUBLIC cxx_std_17)

# Add compile definitions if needed
target_compile_definitions(linearalign PRIVATE
    -DLINEARALIGN_STANDALONE
)
