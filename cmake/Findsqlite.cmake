FIND_PATH(SQLITE_INCLUDE_DIR sqlite3.h 
    PATHS 
    /opt/local/include
    $ENV{HOME}/build/include
    )

FIND_LIBRARY(SQLITE_LIBRARIES NAMES sqlite3 PATHS
    $ENV{HOME}/build/lib
    /opt/local/lib
    /usr/lib64
    )



if (SQLITE_INCLUDE_DIR AND SQLITE_LIBRARIES)
    set (SQLITE_FOUND true)
endif(SQLITE_INCLUDE_DIR AND SQLITE_LIBRARIES)


if (SQLITE_FOUND)
    if (NOT SQLITE_FIND_QUIETLY)
        message(STATUS "Found sqlite: ${SQLITE_LIBRARIES}")
    endif(NOT SQLITE_FIND_QUIETLY)
else (SQLITE_FOUND)
   if (SQLITE_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find sqlite")
    endif(SQLITE_FIND_REQUIRED)
endif(SQLITE_FOUND) 
