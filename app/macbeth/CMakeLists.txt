add_executable(macbeth main.cpp)
    
target_link_libraries(macbeth PUBLIC EXRSpectralImage)

install(TARGETS macbeth RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})