cmake_minimum_required(VERSION 2.6)
set(TARGET SyntheticLightcurves2)
project(${TARGET})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{HOME}/cmake_packages)

find_package(fitsio REQUIRED)
find_package(sqlitepp REQUIRED)
#find_package(Qt4 COMPONENTS QtCore REQUIRED)
find_package(tclap REQUIRED)

#include(${QT_USE_FILE})

include_directories(
    ${${TARGET}_SOURCE_DIR}/include
    ${TCLAP_INCLUDE_DIR}
    ${FITSIO_INCLUDE_DIR}
    ${SQLITEPP_INCLUDE_DIR}
    )

set(SRCS 
    ${${TARGET}_SOURCE_DIR}/src/main.cpp
    )

add_executable(${TARGET} ${SRCS})

target_link_libraries(
    ${TARGET}
    ${FITSIO_LIBRARIES}
    ${SQLITEPP_LIBRARIES}
#    ${QT_LIBRARIES}
    )

