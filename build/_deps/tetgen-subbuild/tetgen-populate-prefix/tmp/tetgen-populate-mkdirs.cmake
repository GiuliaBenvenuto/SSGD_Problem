# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build/_deps/tetgen-src"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build/_deps/tetgen-build"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build/_deps/tetgen-subbuild/tetgen-populate-prefix"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build/_deps/tetgen-subbuild/tetgen-populate-prefix/tmp"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build/_deps/tetgen-subbuild/tetgen-populate-prefix/src/tetgen-populate-stamp"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build/_deps/tetgen-subbuild/tetgen-populate-prefix/src"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build/_deps/tetgen-subbuild/tetgen-populate-prefix/src/tetgen-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build/_deps/tetgen-subbuild/tetgen-populate-prefix/src/tetgen-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build/_deps/tetgen-subbuild/tetgen-populate-prefix/src/tetgen-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
