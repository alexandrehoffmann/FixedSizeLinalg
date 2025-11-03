cmake_minimum_required(VERSION 3.14)
project(LightweightNumericalIntegrationTools VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address")

# List of compiler flags
set(FSLinalg_COMPILE_WARNINGS
    -fopenmp-simd
    -Wundef
    -Wvarargs
    -Wall
    -Wextra
    -Winit-self
    -Wpedantic
    -Werror
    -Wconversion
    -Wuninitialized
    -Wmissing-declarations
    -Wrange-loop-construct
    -Wsign-conversion
    -Wshadow
    -Wnon-virtual-dtor
    -Wold-style-cast
    -Wcast-align
    -Woverloaded-virtual
    -Wnull-dereference
    -Wformat=2
    -flax-vector-conversions
    -pedantic
    -pedantic-errors
	-Wno-error=array-bounds
)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)

    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

# === Options ===
option(FSLinalg_BUILD_DEMO  "Build demo executable" OFF)
option(FSLinalg_BUILD_DOC   "Build Doxygen documentation" OFF)
option(FSLinalg_BUILD_TESTS "Build unit tests" OFF)

# === Dependencies ===
find_package(fmt REQUIRED)

# === Doxygen ===
if(FSLinalg_BUILD_DOC)
    find_package(Doxygen REQUIRED dot)

    # Set default values for Doxyfile
    set(DOXYGEN_PROJECT_NAME "Fixed Size Lightweight Linear algebra library (FSLinalg)")
    set(DOXYGEN_INPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
    set(DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/doc")

    set(DOXYGEN_FILE ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${DOXYGEN_FILE} @ONLY)

    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_FILE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )
endif()

# === Demo subdir ===
if(FSLinalg_BUILD_DEMO)
    add_subdirectory(demo)
endif()

# === Tests ===
if(FSLinalg_BUILD_TESTS)
    enable_testing()

    # Use FetchContent to get GoogleTest if not installed
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    )
    # Prevent GoogleTest from overriding our compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    add_subdirectory(tests)

endif()

# === Library Target ===
add_library(FSLinalg INTERFACE)

target_include_directories(FSLinalg
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(FSLinalg INTERFACE fmt::fmt)

target_compile_options(FSLinalg INTERFACE
    $<$<COMPILE_LANGUAGE:CXX>:${FSLinalg_COMPILE_WARNINGS}>
)

# === Installation ===
include(GNUInstallDirs)

install(DIRECTORY include/FSLinalg DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS FSLinalg
    EXPORT FSLinalgTargets
)

install(EXPORT FSLinalgTargets
    FILE FSLinalgTargets.cmake
    NAMESPACE FSLinalg::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FSLinalg
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/FSLinalgConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FSLinalgConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/FSLinalgConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FSLinalg
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/FSLinalgConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/FSLinalgConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FSLinalg
)

# === pkg-config ===
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FSLinalg.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/FSLinalg.pc
    @ONLY
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/FSLinalg.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)
