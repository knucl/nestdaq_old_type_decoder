cmake_minimum_required(VERSION 3.16)
project(decode)

#if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
#  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" CACHE PATH "Install prefix" FORCE)
#endif()

set(DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(EXECUTABLE_OUTPUT_PATH ${DIR}/install/bin)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT)
if(ROOT_FOUND)
  include(${ROOT_USE_FILE})
endif()

include_directories(
  ${ROOT_INCLUDE_DIRS}
  include
)

link_libraries(
  ${ROOT_LIBRARIES}
)

add_executable(decode       usr/decode.cc)
add_executable(decode_tmp   usr/decode_tmp.cc)
add_executable(decode_slice usr/decode_slice.cc)
add_executable(calc_rate    usr/calc_rate.cc)

add_executable(decode_honda
  usr/decode_honda.cc
  src/NodeUnpacker.cc
  src/LeafProcessor.cc
  src/NodeHeaderDecoder.cc
  src/LeafDecoder.cc
)
target_include_directories(decode_honda PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(decode_honda PRIVATE
  ${ROOT_LIBRARIES}
)

add_executable(check_module_sync   usr/260302/check_module_sync.cc)
add_executable(check_time_interval usr/260303/check_time_interval.cc)
add_executable(check_filter_timing usr/260406/check_filter_timing.cc)
