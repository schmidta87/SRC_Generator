# CMakeLists for simulator library


set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_library(detectors detectors.cpp)

