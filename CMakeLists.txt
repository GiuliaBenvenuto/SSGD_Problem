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
# cmake_minimum_required(VERSION 3.2)
# project(SSGD_dist)

# # GUI application
# add_executable(${PROJECT_NAME} main.cpp solving_ssgd.h solving_ssgd.cpp)
# set_target_properties(${PROJECT_NAME}     PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# # Batch application
# add_executable(batch_app batch.cpp solving_ssgd.h solving_ssgd.cpp) # Assuming batch.cpp includes the same headers and uses the same source files.
# set_target_properties(batch_app PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# # CINOLIB
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

# # Linking libraries common to both applications
# foreach(app IN ITEMS ${PROJECT_NAME} batch_app)
#     target_link_libraries(${app}
#         cinolib
#         graph-based
#         VTP_wrapper
#         Tret
#         geo_central_wrapper
#     )
# endforeach()

# # OpenGL
# find_package(OpenGL REQUIRED)
# include_directories(${OPENGL_INCLUDE_DIR})
# foreach(app IN ITEMS ${PROJECT_NAME} batch_app)
#     target_link_libraries(${app} ${OPENGL_gl_LIBRARY})
# endforeach()

# # Matlab
# set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "Build architectures for Mac OS X" FORCE)
# include_directories("/Applications/MATLAB_R2024a.app/extern/include")
# link_directories("/Applications/MATLAB_R2024a.app/bin/maca64")
# foreach(app IN ITEMS ${PROJECT_NAME} batch_app)
#     target_link_libraries(${app}
#         "/Applications/MATLAB_R2024a.app/extern/bin/maca64/libMatlabDataArray.dylib"
#         "/Applications/MATLAB_R2024a.app/extern/bin/maca64/libMatlabEngine.dylib"
#     )
# endforeach()


cmake_minimum_required(VERSION 3.2)
project(SSGD_dist)

# GUI application
add_executable(${PROJECT_NAME} main.cpp solving_ssgd.h solving_ssgd.cpp)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# SMAPE application
add_executable(smape_app smape.cpp solving_ssgd.h solving_ssgd.cpp)
set_target_properties(smape_app PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# VTP sampling
add_executable(vtp_sampling_app vtp_sampling.cpp solving_ssgd.h solving_ssgd.cpp)
set_target_properties(vtp_sampling_app PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# esp2 sampling
add_executable(esp2_app esp2.cpp solving_ssgd.h solving_ssgd.cpp)
set_target_properties(esp2_app PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# Check Bob Geotangle smape
add_executable(check_bob_app check_bob.cpp solving_ssgd.h solving_ssgd.cpp)
set_target_properties(check_bob_app PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# Smape batch
add_executable(smape_batch_app smape_batch.cpp solving_ssgd.h solving_ssgd.cpp)
set_target_properties(smape_batch_app PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# Write ground truth
add_executable(write_gt_csv write_gt_csv.cpp solving_ssgd.h solving_ssgd.cpp)
set_target_properties(write_gt_csv PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# Save distances to file
add_executable(save_dist_app save_dist.cpp solving_ssgd.h solving_ssgd.cpp)
set_target_properties(save_dist_app PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# Mesh Analysis
add_executable(mesh_analysis_app mesh_analysis.cpp solving_ssgd.h solving_ssgd.cpp)
set_target_properties(mesh_analysis_app PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# CINOLIB
set(cinolib_DIR /Users/giuliabenvenuto/Documents/GitHub/cinolib)
set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)
find_package(cinolib REQUIRED)

# Common include directories and subdirectories
include_directories(
    ${cinolib_DIR}/include
    ${CMAKE_SOURCE_DIR}/libs
    ${CMAKE_SOURCE_DIR}/libs/Eigen
    ${CMAKE_SOURCE_DIR}/libs/geometry-central/include
)
add_subdirectory(SSGD_methods)
add_subdirectory(libs)

# Linking libraries common to all applications
foreach(app IN ITEMS ${PROJECT_NAME} smape_app vtp_sampling_app esp2_app check_bob_app smape_batch_app write_gt_csv save_dist_app mesh_analysis_app)
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
foreach(app IN ITEMS ${PROJECT_NAME} smape_app vtp_sampling_app esp2_app check_bob_app smape_batch_app write_gt_csv save_dist_app mesh_analysis_app)
    target_link_libraries(${app} ${OPENGL_gl_LIBRARY})
endforeach()

# Matlab (assuming the same linking for all applications)
set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "Build architectures for Mac OS X" FORCE)
include_directories("/Applications/MATLAB_R2024a.app/extern/include")
link_directories("/Applications/MATLAB_R2024a.app/bin/maca64")
foreach(app IN ITEMS ${PROJECT_NAME} smape_app vtp_sampling_app esp2_app check_bob_app smape_batch_app write_gt_csv save_dist_app mesh_analysis_app)
    target_link_libraries(${app}
        "/Applications/MATLAB_R2024a.app/extern/bin/maca64/libMatlabDataArray.dylib"
        "/Applications/MATLAB_R2024a.app/extern/bin/maca64/libMatlabEngine.dylib"
    )
endforeach()
