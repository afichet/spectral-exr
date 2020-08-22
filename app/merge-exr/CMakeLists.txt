cmake_minimum_required(VERSION 3.1.1)
project(SpectralImage)

add_executable(merge-exr main.cpp)
    
target_link_libraries(merge-exr PUBLIC EXRSpectralImage)
