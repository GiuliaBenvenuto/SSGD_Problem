# ----- Without batch -----
# cmake_minimum_required(VERSION 3.2)
# project(SSGD_dist)

# add_executable(${PROJECT_NAME} main.cpp solving_ssgd.h solving_ssgd.cpp)
# set_target_properties(${PROJECT_NAME}     PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# set(cinolib_DIR /Users/giuliabenvenuto/Documents/GitHub/cinolib)
# set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)
# find_package(cinolib REQUIRED)

# target_include_directories(${PROJECT_NAME} PUBLIC ${cinolib_DIR}/include)
# target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/libs)
# include_directories(${CMAKE_SOURCE_DIR}/libs/Eigen)

# set(geometry-central_DIR ${CMAKE_SOURCE_DIR}/libs/geometry-central)
# include_directories(${geometry-central_DIR}/include)
# #include_directories(${CMAKE_SOURCE_DIR}/libs/geometry-central/include)
# add_subdirectory(SSGD_methods)
# add_subdirectory(libs)

# target_link_libraries(${PROJECT_NAME} 
#     cinolib 
#     graph-based 
#     VTP_wrapper 
#     Tret 
#     geo_central_wrapper
#     )

# # added for OpenGL
# find_package(OpenGL REQUIRED)
# include_directories(${OPENGL_INCLUDE_DIR})
# target_link_libraries(${PROJECT_NAME} ${OPENGL_gl_LIBRARY})

# # Matlab
# set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "Build architectures for Mac OS X" FORCE)
# include_directories("/Applications/MATLAB_R2024a.app/extern/include")
# link_directories("/Applications/MATLAB_R2024a.app/bin/maca64")
# target_link_libraries(${PROJECT_NAME} 
# "/Applications/MATLAB_R2024a.app/extern/bin/maca64/libMatlabDataArray.dylib"
# "/Applications/MATLAB_R2024a.app/extern/bin/maca64/libMatlabEngine.dylib")


# ----- With batch -----
cmake_minimum_required(VERSION 3.2)
project(SSGD_dist)

# GUI application
add_executable(${PROJECT_NAME} main.cpp solving_ssgd.h solving_ssgd.cpp)
set_target_properties(${PROJECT_NAME}     PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# Batch application
add_executable(batch_app batch.cpp solving_ssgd.h solving_ssgd.cpp) # Assuming batch.cpp includes the same headers and uses the same source files.
set_target_properties(batch_app PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# CINOLIB
set(cinolib_DIR /Users/giuliabenvenuto/Documents/GitHub/cinolib)
set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)
find_package(cinolib REQUIRED)

target_include_directories(${PROJECT_NAME} PUBLIC ${cinolib_DIR}/include)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/libs)
include_directories(${CMAKE_SOURCE_DIR}/libs/Eigen)
set(geometry-central_DIR ${CMAKE_SOURCE_DIR}/libs/geometry-central)
include_directories(${geometry-central_DIR}/include)
#include_directories(${CMAKE_SOURCE_DIR}/libs/geometry-central/include)
add_subdirectory(SSGD_methods)
add_subdirectory(libs)

# Linking libraries common to both applications
foreach(app IN ITEMS ${PROJECT_NAME} batch_app)
    target_link_libraries(${app}
        cinolib
        graph-based
        VTP_wrapper
        Tret
        geo_central_wrapper
    )
endforeach()

# OpenGL
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})
foreach(app IN ITEMS ${PROJECT_NAME} batch_app)
    target_link_libraries(${app} ${OPENGL_gl_LIBRARY})
endforeach()

# Matlab
set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "Build architectures for Mac OS X" FORCE)
include_directories("/Applications/MATLAB_R2024a.app/extern/include")
link_directories("/Applications/MATLAB_R2024a.app/bin/maca64")
foreach(app IN ITEMS ${PROJECT_NAME} batch_app)
    target_link_libraries(${app}
        "/Applications/MATLAB_R2024a.app/extern/bin/maca64/libMatlabDataArray.dylib"
        "/Applications/MATLAB_R2024a.app/extern/bin/maca64/libMatlabEngine.dylib"
    )
endforeach()