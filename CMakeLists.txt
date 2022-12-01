cmake_minimum_required(VERSION 3.2)

project(KuramotoUKF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")

file(GLOB kukf_SRC CONFIGURE_DEPENDS "src/*.h" "src/*.cpp")

add_executable(kukf ${kukf_SRC})


