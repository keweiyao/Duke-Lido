cmake_minimum_required(VERSION 2.8)
project(Lido CXX)

#-----------------#
# Options         #
#-----------------#
option(pythia8 "whether use pythia8 and comiple related examples" on)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" CXX11)
if(NOT CXX11)
  message(FATAL_ERROR "Your compiler does not support C++11")
endif()

# Find and use Boost.
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost 1.60 REQUIRED COMPONENTS filesystem program_options system log log_setup thread REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

# find and use gsl
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
FIND_PACKAGE(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmakemodules")

if(pythia8)
	# Find and use Pythia
	find_package(Pythia8)
	if (${PYTHIA8_FOUND})
	  include_directories(${PYTHIA8_INCLUDE_DIR})
	  Message(STATUS "Pythia8 Include dir : " ${PYTHIA8_INCLUDE_DIR})
	endif()
endif(pythia8)

find_package(FastJet)
if (${FASTJET_FOUND})
  include_directories(${FASTJET_INCLUDE_DIR})
  Message(STATUS "FastJet Include dir : " ${FASTJET_INCLUDE_DIR})
endif()

# Find the HDF5 C++ library.  The CMake FindHDF5 module does not work reliably,
# e.g. sometimes it links unnecessary libraries, other times it fails to include
# libz.  The following works on all systems Jonah Bernhard has tested 
# (Ubuntu, Fedora, Arch, Cray Linux at NERSC, RHEL on OSG).

# Read HDF5_ROOT either as a CMake option or environment variable.
if(NOT HDF5_ROOT)
  set(HDF5_ROOT $ENV{HDF5_ROOT})
endif()

# If HDF5_ROOT was set, search for the header and libraries in that location,
# otherwise try to extract the locations from the h5c++ compiler wrapper.
if(HDF5_ROOT)
  set(HDF5_INCLUDE_HINT "${HDF5_ROOT}/include")
  set(HDF5_LIBRARY_HINT "${HDF5_ROOT}/lib")
else()
  find_program(H5CXX h5c++)
  if(H5CXX)
    if(STATIC)
      set(H5CXX_ARG "-noshlib")
    else()
      set(H5CXX_ARG "-shlib")
    endif()
    # parse output to extract include and library dirs (-I and -L flags)
    execute_process(COMMAND ${H5CXX} -show ${H5CXX_ARG} OUTPUT_VARIABLE H5CXX_OUTPUT)
    foreach(TYPE "INCLUDE" "LIBRARY")
      string(SUBSTRING ${TYPE} 0 1 FLAG)
      set(FLAG " -${FLAG}")
      string(REGEX MATCH "${FLAG}[^ ]+" FLAG_HINT ${H5CXX_OUTPUT})
      string(REPLACE "${FLAG}" "" HDF5_${TYPE}_HINT "${FLAG_HINT}")
    endforeach()
  endif()
endif()

# Find all required HDF5 components.
find_path(HDF5_INCLUDE_DIRS H5Cpp.h HINTS ${HDF5_INCLUDE_HINT})
find_library(HDF5_CXX_LIB hdf5_cpp HINTS ${HDF5_LIBRARY_HINT})
find_library(HDF5_C_LIB hdf5 HINTS ${HDF5_LIBRARY_HINT})
find_library(ZLIB z)

# Group HDF5 components by whether they were found or not.
set(HDF5_VARS_FOUND "")
set(HDF5_VARS_MISSING "")
set(HDF5_FOUND TRUE)
foreach(VAR HDF5_INCLUDE_DIRS HDF5_CXX_LIB HDF5_C_LIB ZLIB)
  if(${VAR})
    set(HDF5_VARS_FOUND "${HDF5_VARS_FOUND} ${${VAR}}")
  else()
    set(HDF5_FOUND FALSE)
    set(HDF5_VARS_MISSING "${HDF5_VARS_MISSING} ${VAR}")
  endif()
endforeach()

# Report HDF5 status and set relevant variables if found.
if(HDF5_FOUND)
  set(HDF5_LIBRARIES ${HDF5_CXX_LIB} ${HDF5_C_LIB} ${ZLIB})
  include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})
  add_definitions(-DTRENTO_HDF5)
  message(STATUS "Found HDF5:${HDF5_VARS_FOUND}")
else()
  message(STATUS "HDF5 not found (missing:${HDF5_VARS_MISSING}, found:${HDF5_VARS_FOUND})")
endif()

set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 1)
set(PROJECT_VERSION_PATCH 1)
set(PROJECT_VERSION
  "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")

# default build type: Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# default install prefix: ~/.local
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local"
      CACHE PATH "Install path prefix, prepended onto install directories."
      FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

#string(APPEND CMAKE_CXX_FLAGS " -march=native")
# By default, optimize for the system's native architecture.  Disable via the
# NATIVE option.  In addition, detect if another architecture flag is already
# set and do not override it.
option(NATIVE "compile for native architecture" ON)
if(NATIVE)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    # Intel compiler: search for -m, -x, -ax flags; if not found add -xHost.
    if(NOT CMAKE_CXX_FLAGS MATCHES "(^| )-(m|x|ax)[^ ]+")
      string(APPEND CMAKE_CXX_FLAGS " -xHost")
    endif()
  else()
    # All other compilers (gcc-like): search for -m flags;
    # if not found add -march=native.
    if(NOT CMAKE_CXX_FLAGS MATCHES "(^| )-m[^ ]+")
      string(APPEND CMAKE_CXX_FLAGS " -march=native")
    endif()
  endif()
endif()

#string(APPEND CMAKE_CXX_FLAGS " -Wall -Wextra")

option(MORE_WARNINGS "enable more compiler warnings" OFF)
if(MORE_WARNINGS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # adapted from http://stackoverflow.com/a/9862800
  string(APPEND CMAKE_CXX_FLAGS " -Wpedantic -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wsign-conversion -Wsign-promo -Wstrict-overflow=2 -Wundef")
endif()

# disable silly clang warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  string(APPEND CMAKE_CXX_FLAGS " -Wno-missing-braces -Wno-c++11-narrowing")
endif()

string(APPEND CMAKE_CXX_FLAGS_DEBUG " -Werror")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  string(APPEND CMAKE_CXX_FLAGS_DEBUG " -Og")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  # Check for sanitizer support.
  # (see e.g. https://code.google.com/p/address-sanitizer/wiki/AddressSanitizer)
  # check_cxx_compiler_flag() doesn't work here because it doesn't set the test
  # flag at link time, so use check_cxx_source_compiles() instead.
  include(CheckCXXSourceCompiles)
  set(CMAKE_REQUIRED_QUIET TRUE)

  # First try both AddressSanitizer and UndefinedBehaviorSanitizer.
  # Disable the vptr sanitizer because it throws errors with boost::any.
  set(CMAKE_REQUIRED_FLAGS "-fsanitize=address,undefined -fno-sanitize=vptr")
  check_cxx_source_compiles("int main() { return 0; }" AddrUndefSanitizer)

  # Lacking that, try AddressSanitizer only.
  if(NOT AddrUndefSanitizer)
    set(CMAKE_REQUIRED_FLAGS "-fsanitize=address")
    check_cxx_source_compiles("int main() { return 0; }" AddrSanitizer)
    if(NOT AddrSanitizer)
      # out of luck...
      set(CMAKE_REQUIRED_FLAGS "")
    endif()
  endif()
endif()

string(APPEND CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_REQUIRED_FLAGS} --coverage")

set(LIBRARY_NAME "lib${PROJECT_NAME}")
add_subdirectory(src)
include_directories(src)

set(LIBRARY_JETFIND "libjetfind")
add_subdirectory(jet_analysis)
include_directories(jet_analysis)


foreach(App "Lido-TabGen" "Lido-X-tests" "Lido-R-tests" "Lido-Incoherent-tests")
	add_executable(${App} ./TestMains/${App}.cpp)
	target_link_libraries(${App} ${LIBRARY_NAME} ${GSL_LIBRARIES} ${GSLCALAS_LIBRARIES} ${HDF5_LIBRARIES} ${Boost_LIBRARIES} -lpthread)
	install(TARGETS ${App} DESTINATION bin)
endforeach()

if(pythia8)
	foreach(App "Lido2DHydro" "Lido_pp")
	add_executable(${App} ./JetMains/${App}.cpp)
	target_link_libraries(${App} ${LIBRARY_NAME} ${LIBRARY_JETFIND} ${PYTHIA8_LIBRARIES} ${FASTJET_LIBRARIES} ${GSL_LIBRARIES} ${GSLCALAS_LIBRARIES} ${HDF5_LIBRARIES} ${Boost_LIBRARIES} -pthread -lpthread -ldl)
	install(TARGETS ${App} DESTINATION bin)
endforeach()

endif(pythia8)

#install(FILES ./settings/lido_settings.xml DESTINATION share)

