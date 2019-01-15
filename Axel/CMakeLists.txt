# CMakeLists for Axel's analysis projects

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_executable(fomin fomin.cpp)
target_link_libraries(fomin ${ROOT_LIBRARIES} )

