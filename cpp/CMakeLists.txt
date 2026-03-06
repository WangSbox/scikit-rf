cmake_minimum_required(VERSION 3.10)
project(skrf_cpp_project VERSION 0.1 LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(${CMAKE_SOURCE_DIR}/cpp/include)

add_executable(skrf_sample src/main.cpp src/network.cpp src/mathFunctions.cpp src/util.cpp src/touchstone.cpp)
target_include_directories(skrf_sample PRIVATE ${CMAKE_SOURCE_DIR}/cpp/include)
target_link_libraries(skrf_sample PRIVATE Eigen3::Eigen)
