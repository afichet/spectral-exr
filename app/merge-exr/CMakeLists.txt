cmake_minimum_required(VERSION 3.1.1)
project(SpectralImage)

# We need C++17 for filesystem library
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(merge-exr main.cpp)
    
target_link_libraries(merge-exr PUBLIC EXRSpectralImage)
