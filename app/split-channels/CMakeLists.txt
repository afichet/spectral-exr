add_executable(split-channels main.cpp)
    
target_link_libraries(split-channels PUBLIC EXRSpectralImage)

install(TARGETS split-channels RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})