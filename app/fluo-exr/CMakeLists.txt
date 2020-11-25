add_executable(fluo-exr main.cpp)
    
target_link_libraries(fluo-exr PUBLIC EXRSpectralImage)

install(TARGETS fluo-exr RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})