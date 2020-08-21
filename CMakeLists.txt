cmake_minimum_required(VERSION 2.8)
project(SpectralImage)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package( OpenEXR REQUIRED)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (OpenEXR_FOUND)
    set(SOURCES
        ${PROJECT_SOURCE_DIR}/src/SpectralImage.cpp
        ${PROJECT_SOURCE_DIR}/src/EXRSpectralImage.cpp

        # Optional bi spectral variants
        ${PROJECT_SOURCE_DIR}/src/BiSpectralImage.cpp
        ${PROJECT_SOURCE_DIR}/src/EXRBiSpectralImage.cpp

        ${PROJECT_SOURCE_DIR}/src/main.cpp
    )

    add_executable(
        SpectralImage
        ${SOURCES}
    )

    target_include_directories(
        SpectralImage PUBLIC
        ${PROJECT_SOURCE_DIR}/src
        ${OpenEXR_INCLUDE_DIR}
    )
    
    target_link_libraries(
        SpectralImage PUBLIC
        ${OpenEXR_LIBRARIES}
    )
else()
    message("You need OpenEXR libraries to compile this code.")
endif()