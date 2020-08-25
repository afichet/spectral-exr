cmake_minimum_required(VERSION 3.1.1)
project(SpectralImage)

add_executable(fluo-exr main.cpp)
    
target_link_libraries(fluo-exr PUBLIC EXRSpectralImage)
