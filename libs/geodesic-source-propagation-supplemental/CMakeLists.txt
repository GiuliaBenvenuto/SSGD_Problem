cmake_minimum_required(VERSION 3.0)

project(InstantGeodesics)
add_executable(${PROJECT_NAME} "main.cc")

set_target_properties(${PROJECT_NAME}     PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_STANDARD 17)

if (MSVC)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE
        -ffast-math
        -march=native
        )
endif()
