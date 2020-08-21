cmake_minimum_required(VERSION 3.1.1)
project(SpectralImage)

add_executable(split-channels main.cpp)
    
target_link_libraries(split-channels PUBLIC EXRSpectralImage)
