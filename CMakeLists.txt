cmake_minimum_required(VERSION 3.16)

project(lfsr128sum LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options("-msse4.1")
add_executable(lfsr128sum main.cpp
    lfsr_hash.h tests.h)

include(GNUInstallDirs)
install(TARGETS lfsr128sum
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
