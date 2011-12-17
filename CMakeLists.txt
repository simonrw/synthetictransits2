cmake_minimum_required(VERSION 2.6)
set(TARGET SyntheticTransits)
project(${TARGET})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${${TARGET}_SOURCE_DIR}/cmake)

find_package(fitsio REQUIRED)
find_package(sqlitepp REQUIRED)
#find_package(Qt4 COMPONENTS QtCore REQUIRED)
find_package(tclap REQUIRED)

#include(${QT_USE_FILE})


add_subdirectory(Modelgen)

include_directories(
    ${${TARGET}_SOURCE_DIR}/include
    ${MODELGEN_INCLUDE_DIR}
    ${TCLAP_INCLUDE_DIR}
    ${FITSIO_INCLUDE_DIR}
    ${SQLITEPP_INCLUDE_DIR}
    )

set(SRCS 
    ${${TARGET}_SOURCE_DIR}/src/main.cpp
    ${${TARGET}_SOURCE_DIR}/src/timer.C
    ${${TARGET}_SOURCE_DIR}/src/FitsObject.cpp
    )

set(HDRS
    ${${TARGET}_SOURCE_DIR}/include/timer.h
    ${${TARGET}_SOURCE_DIR}/include/FitsObject.h
    ${${TARGET}_SOURCE_DIR}/include/ObjectSkipDefs.h
    )

add_executable(${TARGET} ${SRCS} ${HDRS})
add_dependencies(${TARGET} ${MODELTARGET})

target_link_libraries(
    ${TARGET}
    ${FITSIO_LIBRARIES}
    ${SQLITEPP_LIBRARIES}
    ${MODELGEN_LIBRARIES}
#    ${QT_LIBRARIES}
    )

install(TARGETS ${TARGET} DESTINATION bin)
set_property(TARGET ${TARGET} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
