# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/indirect_predicates-src"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/indirect_predicates-build"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/indirect_predicates-subbuild/indirect_predicates-populate-prefix"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/indirect_predicates-subbuild/indirect_predicates-populate-prefix/tmp"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/indirect_predicates-subbuild/indirect_predicates-populate-prefix/src/indirect_predicates-populate-stamp"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/indirect_predicates-subbuild/indirect_predicates-populate-prefix/src"
  "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/indirect_predicates-subbuild/indirect_predicates-populate-prefix/src/indirect_predicates-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/indirect_predicates-subbuild/indirect_predicates-populate-prefix/src/indirect_predicates-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/indirect_predicates-subbuild/indirect_predicates-populate-prefix/src/indirect_predicates-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
