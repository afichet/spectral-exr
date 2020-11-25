add_executable(export-spectrum main.cpp)
    
target_link_libraries(export-spectrum PUBLIC EXRSpectralImage)

install(TARGETS export-spectrum RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})