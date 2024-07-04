add_library(geo_central_wrapper gc_wrapper.h gc_wrapper.cpp)
set_target_properties(geo_central_wrapper PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

#set(geometry-central_DIR ${CMAKE_SOURCE_DIR}/geometry-central)
#set(geometry-central_DIR ${CMAKE_SOURCE_DIR}/../../libs/geometry-central)
#set(geometry-central_DIR ${CMAKE_SOURCE_DIR}/../libs/geometry-central)




target_include_directories(geo_central_wrapper PUBLIC ${cinolib_DIR}/include)
target_include_directories(geo_central_wrapper PUBLIC ${geometry-central_DIR}/include)

target_link_libraries(geo_central_wrapper cinolib geometry-central)
