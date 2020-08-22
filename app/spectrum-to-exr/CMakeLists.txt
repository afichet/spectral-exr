cmake_minimum_required(VERSION 3.1.1)
project(SpectralImage)

add_executable(spectrum-to-exr main.cpp)
    
target_link_libraries(spectrum-to-exr PUBLIC EXRSpectralImage)
