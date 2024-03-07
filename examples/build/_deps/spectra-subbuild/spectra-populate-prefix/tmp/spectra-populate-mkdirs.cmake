# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/spectra-src"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/spectra-build"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/spectra-subbuild/spectra-populate-prefix"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/spectra-subbuild/spectra-populate-prefix/tmp"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/spectra-subbuild/spectra-populate-prefix/src/spectra-populate-stamp"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/spectra-subbuild/spectra-populate-prefix/src"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/spectra-subbuild/spectra-populate-prefix/src/spectra-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/spectra-subbuild/spectra-populate-prefix/src/spectra-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/spectra-subbuild/spectra-populate-prefix/src/spectra-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
