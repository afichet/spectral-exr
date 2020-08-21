cmake_minimum_required(VERSION 3.1.1)
project(SpectralImage)

add_executable(SpectralImage
    main.cpp
)
    
target_link_libraries(SpectralImage PUBLIC EXRSpectralImage)
