cmake_minimum_required(VERSION 3.19)

project(VBMicrolensing LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(MINIMUM_GCC_VERSION "9.0")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS MINIMUM_GCC_VERSION)
        message(FATAL_ERROR "GCC must be ${MINIMUM_GCC_VERSION} or greater.")
    endif()
endif()

# Find the module development requirements (requires FindPython from 3.17 or
# scikit-build-core's built-in backport)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Add a library using FindPython's tooling (pybind11 also provides a helper like
# this)
python_add_library(VBMicrolensing MODULE VBMicrolensing/lib/python_bindings.cpp VBMicrolensing/lib/VBMicrolensingLibrary.cpp VBMicrolensing/lib/VBMicrolensingLibrary.h)
target_link_libraries(VBMicrolensing PRIVATE pybind11::headers)

# The install directory is the output (wheel) directory
install(TARGETS VBMicrolensing DESTINATION VBMicrolensing)
