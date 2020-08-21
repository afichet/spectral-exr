cmake_minimum_required(VERSION 3.1.1)
project(SpectralImage)

#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#find_package( OpenEXR REQUIRED)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_subdirectory(lib)
add_subdirectory(app)
