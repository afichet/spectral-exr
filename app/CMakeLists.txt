cmake_minimum_required(VERSION 3.1.1)
project(SpectralImage)

add_subdirectory(macbeth)
add_subdirectory(spectrum-to-exr)

add_subdirectory(split-channels)
add_subdirectory(merge-exr)

