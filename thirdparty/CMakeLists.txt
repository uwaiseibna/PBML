# Thirdparty libraries
include(FetchContent)
include(ExternalProject)

# Log directory for external content builds
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/logs)
set(INSTALL_FILE ${CMAKE_CURRENT_BINARY_DIR}/install.cmake)

# Set environment
list(APPEND CMAKE_LIBRARY_PATH "${CMAKE_CURRENT_BINARY_DIR}/lib")
list(APPEND CMAKE_INCLUDE_PATH "${CMAKE_CURRENT_BINARY_DIR}/include")
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_BINARY_DIR}")

# External Project resources
# ------------------------------------------------------------------------------

# SDSL (includes divsufsort and divsufsort64)
find_library(SDSL_LIB sdsl HINTS ${CMAKE_CURRENT_BINARY_DIR}/lib PATHS ${CONDA_PREFIX}/lib)
find_path(SDSL_SRC sdsl HINTS ${CMAKE_CURRENT_BINARY_DIR}/include PATHS ${CONDA_PREFIX}/include)

if(NOT SDSL_LIB)
    message(STATUS "sdsl library not found. Building as an external content")
    FetchContent_Declare(
        sdsl
        GIT_REPOSITORY "https://github.com/simongog/sdsl-lite.git"
        GIT_TAG "v2.1.1"
    )

    FetchContent_MakeAvailable(sdsl)
    
    # SDSL's CMakeLists.txt already creates targets: sdsl, divsufsort, divsufsort64
    # No need to create IMPORTED targets - they already exist from the subdirectory
    
else()
    message(STATUS "sdsl library found at ${SDSL_LIB}.")
    message(STATUS "sdsl sources found at ${SDSL_SRC}.")
    
    # Find divsufsort libraries (required when using system SDSL)
    find_library(DIVSUFSORT_LIB divsufsort HINTS ${CMAKE_CURRENT_BINARY_DIR}/lib PATHS ${CONDA_PREFIX}/lib REQUIRED)
    find_library(DIVSUFSORT64_LIB divsufsort64 HINTS ${CMAKE_CURRENT_BINARY_DIR}/lib PATHS ${CONDA_PREFIX}/lib REQUIRED)
    
    # Only create IMPORTED targets if using pre-installed libraries
    add_library(sdsl STATIC IMPORTED GLOBAL)
    set_property(TARGET sdsl PROPERTY IMPORTED_LOCATION ${SDSL_LIB})
    target_include_directories(sdsl INTERFACE ${SDSL_SRC})

    add_library(divsufsort STATIC IMPORTED GLOBAL)
    set_property(TARGET divsufsort PROPERTY IMPORTED_LOCATION ${DIVSUFSORT_LIB})
    
    add_library(divsufsort64 STATIC IMPORTED GLOBAL)
    set_property(TARGET divsufsort64 PROPERTY IMPORTED_LOCATION ${DIVSUFSORT64_LIB})
endif()

## HTSLIB
FetchContent_Declare(
    htslib
    GIT_REPOSITORY https://github.com/samtools/htslib
    GIT_TAG 1.20
)

FetchContent_GetProperties(htslib)
if(NOT htslib_POPULATED)
    FetchContent_Populate(htslib)
    
    set(disable_flags --disable-gcs --disable-s3 --disable-plugins --disable-lzma 
        --disable-bz2 --disable-libcurl --without-libdeflate)

    message(STATUS "HTSLIB_SRC: ${htslib_SOURCE_DIR}")
    message(STATUS "HTSLIB_BIN: ${htslib_BINARY_DIR}")
    
    execute_process(
        COMMAND autoreconf -i
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${htslib_SOURCE_DIR}
        OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib_autoreconf.log
        ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib_autoreconf.log
    )
    if(result)
        message(FATAL_ERROR "autoreconf step for htslib failed: ${result}")
    endif()
    
    execute_process(
        COMMAND ./configure --prefix=${CMAKE_CURRENT_BINARY_DIR} ${disable_flags}
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${htslib_SOURCE_DIR}
        OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib_configure.log
        ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib_configure.log
    )
    if(result)
        message(FATAL_ERROR "configure step for htslib failed: ${result}")
    endif()
    
    execute_process(
        COMMAND make lib-static -j${CMAKE_BUILD_PARALLEL_LEVEL}
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${htslib_SOURCE_DIR}
        OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib_make.log
        ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib_make.log
    )
    if(result)
        message(FATAL_ERROR "make step for htslib failed: ${result}")
    endif()
    
    execute_process(
        COMMAND make install prefix=${CMAKE_CURRENT_BINARY_DIR}
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${htslib_SOURCE_DIR}
        OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib_install.log
        ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib_install.log
    )
    if(result)
        message(FATAL_ERROR "install step for htslib failed: ${result}")
    endif()
    
    set(htslib_SRC ${CMAKE_CURRENT_BINARY_DIR}/include)
    set(htslib_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libhts.a)
endif()

add_library(htslib STATIC IMPORTED GLOBAL)
set_property(TARGET htslib PROPERTY IMPORTED_LOCATION ${htslib_LIB})
target_include_directories(htslib INTERFACE ${htslib_SRC})