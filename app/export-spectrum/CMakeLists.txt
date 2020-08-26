cmake_minimum_required(VERSION 3.1.1)
project(SpectralImage)

add_executable(export-spectrum main.cpp)
    
target_link_libraries(export-spectrum PUBLIC EXRSpectralImage)
