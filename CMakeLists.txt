cmake_minimum_required(VERSION 3.17)
project(heat)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -O3 -mfma -mavx2 -lpython3.8")

include_directories(.)
include_directories(src)
include_directories(/usr/include/python3.8/)
include_directories(/home/alexshokolarov/Documents/Resources/cpplibs/matplotlib-cpp/)

find_package(Python3 COMPONENTS Development)

add_executable(heat
        src/linalg.h
        src/heat.h
        main.cpp)

target_include_directories(heat PRIVATE ${PYTHON3.8_INCLUDE_DIRS})
target_link_libraries(heat Python3::Python)
