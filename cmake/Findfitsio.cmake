FIND_PACKAGE(Threads REQUIRED)

FIND_PATH(FITSIO_INCLUDE_DIR fitsio.h PATHS 
    $ENV{HOME}/build/include
    /opt/local/include
    /usr/local/cfitsio-thread # for wasphead
    )

FIND_LIBRARY(FITSIO_ONLY_LIBRARIES NAMES cfitsio PATHS
    $ENV{HOME}/build/lib
    /opt/local/lib
    /usr/local/cfitsio-thread # for wasphead
    #/sw/lib
    )

if (FITSIO_INCLUDE_DIR AND FITSIO_ONLY_LIBRARIES AND CMAKE_THREAD_LIBS_INIT)
    set (FITSIO_FOUND true)
endif (FITSIO_INCLUDE_DIR AND FITSIO_ONLY_LIBRARIES AND CMAKE_THREAD_LIBS_INIT)

if (FITSIO_FOUND)
    set(FITSIO_LIBRARIES ${FITSIO_ONLY_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    if (NOT FITSIO_FIND_QUIETLY)
        message(STATUS "Found fitsio: ${FITSIO_LIBRARIES}")
    endif(NOT FITSIO_FIND_QUIETLY)
else (FITSIO_FOUND)
   if (FITSIO_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find fitsio")
    endif(FITSIO_FIND_REQUIRED)
endif(FITSIO_FOUND) 
