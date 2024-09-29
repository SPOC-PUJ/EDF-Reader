cmake_minimum_required(VERSION 3.22)
project(DSP)
#set(CMAKE_BUILD_TYPE Debug)
add_executable(EDFReader main.cpp edf.cpp signaldata.cpp)

# Add compiler options for more warnings and checks
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(EDFReader PRIVATE
        -Wall            # Enable all the most common warnings
        -Wextra          # Enable extra warnings
        -Wpedantic       # Enforce strict ISO C++ standards compliance
        -Wshadow         # Warn if a variable declaration shadows one from a parent scope
        -Wnon-virtual-dtor # Warn when classes with virtual functions don't have a virtual destructor
        -Wcast-align     # Warn about potential performance problem casts
        -Woverloaded-virtual # Warn when a function declaration hides a virtual function from a base class
        -Wconversion     # Warn for implicit type conversions that may change a value
        -Wduplicated-cond # Warn if 'if' and 'else' have identical conditions
        -Wnull-dereference # Warn if a null pointer dereference is detected
    )
endif()

add_custom_target(ValgrindCheck
    COMMAND valgrind --leak-check=full --track-origins=yes --show-reachable=yes --error-exitcode=1 ./EDFReader ${CMAKE_CURRENT_BINARY_DIR}/testfiles/test.edf
    DEPENDS EDFReader
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running Valgrind to check for memory leaks"
)


