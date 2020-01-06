# CMakeLists for simulator library

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

# Igor's Hall A Simulator
add_executable(Korover_HallA_sim Korover_HallA_sim.cpp)
target_link_libraries(Korover_HallA_sim ${ROOT_LIBRARIES} detectors)

# Ran's Hall A Simulator
add_executable(Shneor_HallA_sim Shneor_HallA_sim.cpp)
target_link_libraries(Shneor_HallA_sim ${ROOT_LIBRARIES} detectors)

add_library(detectors detectors.cpp)
target_link_libraries(detectors ${ROOT_LIBRARIES} helpers)

