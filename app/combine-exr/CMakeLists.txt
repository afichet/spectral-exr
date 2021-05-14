add_executable(combine-exr main.cpp)
    
target_link_libraries(combine-exr PUBLIC EXRSpectralImage)

install(TARGETS combine-exr RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})