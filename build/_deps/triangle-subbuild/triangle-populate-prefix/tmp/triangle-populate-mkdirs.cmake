# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build/_deps/triangle-src"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build/_deps/triangle-build"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build/_deps/triangle-subbuild/triangle-populate-prefix"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build/_deps/triangle-subbuild/triangle-populate-prefix/tmp"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build/_deps/triangle-subbuild/triangle-populate-prefix/src/triangle-populate-stamp"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build/_deps/triangle-subbuild/triangle-populate-prefix/src"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build/_deps/triangle-subbuild/triangle-populate-prefix/src/triangle-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build/_deps/triangle-subbuild/triangle-populate-prefix/src/triangle-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build/_deps/triangle-subbuild/triangle-populate-prefix/src/triangle-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
