FIND_PATH(SQLITE_INCLUDE_DIR sqlite3.h 
    PATHS 
    /opt/local/include
    $ENV{HOME}/build/include
    )

FIND_LIBRARY(SQLITE3_LIBRARIES NAMES sqlite3 PATHS
    $ENV{HOME}/build/lib
    /opt/local/lib
    NO_DEFAULT_PATH
    )

find_library(SQLITEOLD_LIBRARIES NAMES sqlite PATHS
    $ENV{HOME}/build/lib
    /opt/local/lib
    NO_DEFAULT_PATH
    )


FIND_LIBRARY(DLLIBRARY NAMES dl)

if (SQLITE3_LIBRARIES)
    set(SQLITE_LIBRARIES ${DLLIBRARY} ${SQLITE3_LIBRARIES})
else()
    set(SQLITE_LIBRARIES ${DLLIBRARY} ${SQLITEOLD_LIBRARIES})
endif()



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
