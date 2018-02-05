# - Basic CMake setup
# Check version meets ou requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.3)
project(BiPo)

# Modules use Falaise, so we need to locate this or fail
# Locating Falaise will automatically locate all of its
# dependencies such as Bayeux, ROOT and Boost.
find_package(Falaise REQUIRED)

# Build a dynamic library from our sources
add_library(BiPo SHARED
  BiPo.h BiPo.cc
  )

# Link it to the FalaiseModule library
# This ensures the correct compiler flags, include paths
# and linker flags are applied to our dynamic library.
target_link_libraries(BiPo PUBLIC Falaise::FalaiseModule)

#include_directories(/Users/minotti/SuperNEMO/sw/Falaise/install/include/falaise/snemo)
include_directories($ENV{FALAISE_PATH}/include/falaise/snemo)

