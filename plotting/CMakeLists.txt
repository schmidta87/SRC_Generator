# CMakeLists for simulator library

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

# Igor's Hall A Simulator
add_executable(Korover_HallA_hists Korover_HallA_hists.cpp)
target_link_libraries(Korover_HallA_hists ${ROOT_LIBRARIES} helpers)

# Ran's Hall A Simulator
add_executable(Shneor_HallA_hists Shneor_HallA_hists.cpp)
target_link_libraries(Shneor_HallA_hists ${ROOT_LIBRARIES} helpers)
