add_executable(spectrum-to-exr main.cpp)
    
target_link_libraries(spectrum-to-exr PUBLIC EXRSpectralImage)

install(TARGETS spectrum-to-exr RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})