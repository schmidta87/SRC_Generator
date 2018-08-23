# CMakeLists for SRC_Generator

# Set up
cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)

project (SRC_Generator)

add_executable(gen_weight gen_weight.cpp )
target_link_libraries(gen_weight ${ROOT_LIBRARIES} )
