cmake_minimum_required(VERSION 3.24.2)
project(Tesi)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(Tesi main.cpp
        metriche.cpp metriche.hpp
        fraction.cpp fraction.hpp
        int_mod.cpp int_mod.hpp)
