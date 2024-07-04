add_library(graph-based shortest_path.h shortest_path.cpp geodesic_solver.h geodesic_solver.cpp extended_solver.h extended_solver.cpp)

set_target_properties(graph-based PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

target_include_directories(graph-based PUBLIC ${cinolib_DIR}/include)
target_link_libraries(graph-based cinolib)