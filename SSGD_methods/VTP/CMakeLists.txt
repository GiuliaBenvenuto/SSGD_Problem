add_library(VTP_wrapper vtp_wrapper.h vtp_wrapper.cpp)

set_target_properties(VTP_wrapper PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

target_include_directories(VTP_wrapper PUBLIC ${cinolib_DIR}/include)
target_link_libraries(VTP_wrapper cinolib)