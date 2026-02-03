cmake_minimum_required(VERSION 3.15...3.27)
project(zna_accel LANGUAGES CXX)

# Find Python and nanobind
find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)

# Try to import nanobind from pip install
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE nanobind_ROOT
    RESULT_VARIABLE nanobind_FOUND
)

if(nanobind_FOUND EQUAL 0 AND nanobind_ROOT)
    list(APPEND CMAKE_PREFIX_PATH "${nanobind_ROOT}")
endif()

find_package(nanobind CONFIG REQUIRED)

# Build the extension module
nanobind_add_module(
    _accel
    STABLE_ABI
    NB_STATIC
    src/zna/_accel.cpp
)

# Set C++ standard
target_compile_features(_accel PRIVATE cxx_std_17)

# Optimization flags
if(NOT MSVC)
    target_compile_options(_accel PRIVATE -O3)
else()
    target_compile_options(_accel PRIVATE /O2)
endif()

# Install the module into the zna package
install(TARGETS _accel LIBRARY DESTINATION zna)
