add_library(Tret trettner.h trettner.cpp)

set_target_properties(Tret PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

target_include_directories(Tret PUBLIC ${cinolib_DIR}/include)
target_link_libraries(Tret cinolib)