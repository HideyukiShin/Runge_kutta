cmake_minimum_required(VERSION 3.0)
project(test)
include_directories("include")
find_package(matplotlib_cpp REQUIRED)
find_package(Eigen3 REQUIRED)

add_executable(
  test 
  test.cpp
)

target_include_directories(
  test PRIVATE
  ${matplotlib_cpp_INCLUDE_DIRS}
  ${Python3_INCLUDE_DIRS}
  ${Python3_NumPy_INCLUDE_DIRS}
  ${Eigen3_INCLUDE_DIRS}
)
target_link_libraries(test "${Python3_LIBRARIES}" Eigen3::Eigen)
