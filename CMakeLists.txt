cmake_minimum_required(VERSION 3.16)
project(PBML)
set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)

# Add OpenMP compile flags (THIS WAS MISSING!)
add_compile_options(-Wall -Wextra -O3 -DNDEBUG -march=native ${OpenMP_CXX_FLAGS})

# Configure thirdparty dependencies
set(CMAKE_INSTALL_INCLUDEDIR "include")
add_subdirectory(thirdparty)

# Build executable
add_executable(pbml src/pbml.cpp)

# Link libraries including OpenMP
target_link_libraries(pbml 
    sdsl 
    divsufsort 
    divsufsort64 
    htslib 
    z 
    bz2
    lzma
    curl
    pthread
    m
    OpenMP::OpenMP_CXX  # Better way to link OpenMP
)

# Installation
set(tool_dest "bin")
install(TARGETS pbml DESTINATION "${tool_dest}")