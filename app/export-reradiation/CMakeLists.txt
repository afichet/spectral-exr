add_executable(export-reradiation main.cpp)
    
target_link_libraries(export-reradiation PUBLIC EXRSpectralImage)

install(TARGETS export-reradiation RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})