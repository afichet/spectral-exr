cmake_minimum_required(VERSION 3.1.1)

project(EXRSpectralImage
    VERSION 1.0.0 
    DESCRIPTION "Library for loading and writing OpenEXR spectral images"
    )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# OpenEXR 3.0 provides its own find_package
find_package(OpenEXR CONFIG)
find_package(Imath CONFIG)

if (OpenEXR_FOUND)
    set(OpenEXR_LIBRARIES OpenEXR::OpenEXR Imath::Imath)
else()
    # For OpenEXR 2.0 we rely on our own find_package
    set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
    find_package(OpenEXR REQUIRED)
endif()

if (OpenEXR_FOUND)
    set(PUBLIC_HEADERS
        include/SpectrumType.h
        include/SpectrumAttribute.h
        
        include/SpectralImage.h
        include/EXRSpectralImage.h

        # Optional bi spectral variants
        include/BiSpectralImage.h
        include/EXRBiSpectralImage.h
    )

    add_library(EXRSpectralImage SHARED
        SpectralImage.cpp
        EXRSpectralImage.cpp

        SpectrumConverter.cpp
        SpectrumAttribute.cpp

        # Optional bi spectral variants
        BiSpectralImage.cpp
        EXRBiSpectralImage.cpp
    )

    set_target_properties(EXRSpectralImage PROPERTIES VERSION ${PROJECT_VERSION})
    set_target_properties(EXRSpectralImage PROPERTIES SOVERSION 1)
    set_target_properties(EXRSpectralImage PROPERTIES PUBLIC_HEADER "${PUBLIC_HEADERS}")

    target_include_directories(EXRSpectralImage PUBLIC include)
    target_include_directories(EXRSpectralImage PUBLIC ${OpenEXR_INCLUDE_DIR})
    
    target_link_libraries(EXRSpectralImage PUBLIC ${OpenEXR_LIBRARIES})

    if (MSVC)
        target_compile_options(EXRSpectralImage PUBLIC /W3)
    else()
        target_compile_options(EXRSpectralImage PUBLIC -Wall -Wextra -Wpedantic)
    endif()

    configure_file(EXRSpectralImage.pc.in ${CMAKE_BINARY_DIR}/EXRSpectralImage.pc @ONLY)

    install(TARGETS EXRSpectralImage
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    install(FILES ${CMAKE_BINARY_DIR}/EXRSpectralImage.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
else()
    message("You need OpenEXR libraries to compile this code.")
endif()
