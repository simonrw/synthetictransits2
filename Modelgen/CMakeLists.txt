cmake_minimum_required(VERSION 2.6)
set(TARGET Modelgen)
project(${TARGET})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

include_directories(
    ${${TARGET}_SOURCE_DIR}/include
    )

add_library(${TARGET}
    ${${TARGET}_SOURCE_DIR}/src/GenerateModel.cpp
    )

