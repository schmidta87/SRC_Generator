# CMakeLists for SRC_Generator

# Set up
cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)

project (SRC_Generator)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES} ${input_csv})

add_subdirectory(helpers)
add_subdirectory(simulator)
add_subdirectory(plotting)
add_subdirectory(Andrew)
add_subdirectory(Axel)
add_subdirectory(Mean_Field_Generator)

add_executable(gen_weight gen_weight.cpp Nuclear_Info.cpp Cross_Sections.cpp)
target_link_libraries(gen_weight ${ROOT_LIBRARIES} helpers)

add_executable(smearTree smearTree.cpp)
target_link_libraries(smearTree ${ROOT_LIBRARIES})
