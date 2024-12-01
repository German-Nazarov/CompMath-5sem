cmake_minimum_required(VERSION 3.5.1)
project(Lab-RK VERSION 0.1.0 LANGUAGES C CXX) 

include(CTest)
enable_testing()

add_executable(Lab-RK lab-1-rk.cpp)

set_property(TARGET Lab-RK PROPERTY CXX_STANDARD 17)