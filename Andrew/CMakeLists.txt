# CMakeLists for SRC_Generator


set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_executable(src_reader src_reader.cpp)
target_link_libraries(src_reader ${ROOT_LIBRARIES} )

add_executable(pMissReader pMissReader.cpp)
target_link_libraries(pMissReader ${ROOT_LIBRARIES} )
