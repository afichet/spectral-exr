cmake_minimum_required(VERSION 3.1.1)

project(EXRSpectralImage
    VERSION 1.0.0 
    DESCRIPTION "Library for loading and writing OpenEXR spectral images"
    )

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenEXR REQUIRED)

if (OpenEXR_FOUND)
    set(PUBLIC_HEADERS
        include/SpectrumAttribute.h
        include/SpectralImage.h
        include/EXRSpectralImage.h

        # Optional bi spectral variants
        #include/BiSpectralImage.h
        #include/EXRBiSpectralImage.h
    )

    add_library(EXRSpectralImage SHARED
        SpectralImage.cpp
        EXRSpectralImage.cpp

        SpectrumConverter.cpp
        SpectrumAttribute.cpp

        # Optional bi spectral variants
        #BiSpectralImage.cpp
        #EXRBiSpectralImage.cpp
    )

    set_target_properties(EXRSpectralImage PROPERTIES VERSION ${PROJECT_VERSION})
    set_target_properties(EXRSpectralImage PROPERTIES SOVERSION 1)
    set_target_properties(EXRSpectralImage PROPERTIES PUBLIC_HEADER "${PUBLIC_HEADERS}")

    target_include_directories(EXRSpectralImage PUBLIC include)
    target_include_directories(EXRSpectralImage PUBLIC ${OpenEXR_INCLUDE_DIR})
    
    target_link_libraries(EXRSpectralImage PUBLIC ${OpenEXR_LIBRARIES})

    if(MSVC)
        target_compile_options(EXRSpectralImage PUBLIC /W3)
    else()
        target_compile_options(EXRSpectralImage PUBLIC -Wall -Wextra -Wpedantic)
    endif()

    configure_file(EXRSpectralImage.pc.in EXRSpectralImage.pc @ONLY)

    install(TARGETS EXRSpectralImage
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    install(FILES ${CMAKE_BINARY_DIR}/EXRSpectralImage.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
else()
    message("You need OpenEXR libraries to compile this code.")
endif()