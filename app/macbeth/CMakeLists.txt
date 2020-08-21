cmake_minimum_required(VERSION 3.1.1)
project(SpectralImage)

add_executable(macbeth main.cpp)
    
target_link_libraries(macbeth PUBLIC EXRSpectralImage)
