cmake_minimum_required(VERSION 2.6)
set(TARGET SyntheticLightcurves2)
project(${TARGET})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

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
    ${${TARGET}_SOURCE_DIR}/src/timer.C
    ${${TARGET}_SOURCE_DIR}/src/GenerateModel.cpp
    ${${TARGET}_SOURCE_DIR}/src/FitsObject.cpp
    )

set(HDRS
    ${${TARGET}_SOURCE_DIR}/include/timer.h
    ${${TARGET}_SOURCE_DIR}/include/constants.h
    ${${TARGET}_SOURCE_DIR}/include/Model.h
    ${${TARGET}_SOURCE_DIR}/include/GenerateModel.h
    ${${TARGET}_SOURCE_DIR}/include/FitsObject.h
    ${${TARGET}_SOURCE_DIR}/include/ObjectSkipDefs.h
    )

add_executable(${TARGET} ${SRCS} ${HDRS})

target_link_libraries(
    ${TARGET}
    ${FITSIO_LIBRARIES}
    ${SQLITEPP_LIBRARIES}
#    ${QT_LIBRARIES}
    )

