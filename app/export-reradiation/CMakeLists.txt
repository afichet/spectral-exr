cmake_minimum_required(VERSION 3.1.1)
project(SpectralImage)

add_executable(export-reradiation main.cpp)
    
target_link_libraries(export-reradiation PUBLIC EXRSpectralImage)
