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

# Build the extension module.
#
# No STABLE_ABI here.  It was requested before, but silently did nothing:
# nanobind can only build against the limited API when CMake finds Python's
# Development.SABIModule component, and the find_package above asks for
# Development.Module.  The built artefact was always version-tagged
# (_accel.cpython-3XX-*.so), never abi3.
#
# That is not an oversight to correct — the codec cannot be built on the limited
# API.  It has used PyTuple_SET_ITEM since before the 0.3.4 rewrite, and the
# rewrite also needs PyUnicode_New, PyUnicode_DATA and PyList_SET_ITEM to keep
# from copying every base twice.  Asking for a stable ABI here would simply fail
# to compile, so the request is dropped rather than left misleading.
nanobind_add_module(
    _accel
    NB_STATIC
    src/zna/_accel.cpp
)

# On macOS, nanobind does not link against libpython.  It instead passes the
# linker a curated list of CPython symbols that are allowed to be undefined
# (nanobind/cmake/darwin-ld-cpython.sym), and PyUnicode_New is not on it, so the
# decoder fails to link without this.
#
# The decoder writes bases directly into each str's own buffer; the alternative
# that links unaided (decode to scratch, then PyUnicode_FromStringAndSize) scans
# its input as UTF-8 and copies, costing three passes over every sequence
# instead of one -- measured at a third of the decode win.
#
# Linux and Windows resolve PyUnicode_New normally and need no flag.  If a
# future toolchain rejects -U, this fails as a link error at build time, which
# is the failure mode to want: loud, and impossible to ship past.
if(APPLE)
    target_link_options(_accel PRIVATE "-Wl,-U,_PyUnicode_New")
endif()

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
