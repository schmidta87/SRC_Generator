# CMakeLists for SRC_Generator

# Set up
cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)

project (SRC_Generator)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_subdirectory(Andrew)

add_executable(gen_weight gen_weight.cpp Nuclear_Info.cpp)
target_link_libraries(gen_weight ${ROOT_LIBRARIES} )
