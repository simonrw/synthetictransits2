FIND_PATH(TCLAP_INCLUDE_DIR CmdLine.h PATHS 
    $ENV{HOME}/build/include
    /sw/include
    /opt/local/include
    PATH_SUFFIXES
    tclap
    )


if (TCLAP_INCLUDE_DIR)
    set (TCLAP_FOUND true)
endif(TCLAP_INCLUDE_DIR)

if (TCLAP_FOUND)
    if (NOT TCLAP_FIND_QUIETLY)
        message(STATUS "Found TCLAP: ${TCLAP_INCLUDE_DIR}")
    endif(NOT TCLAP_FIND_QUIETLY)
else (TCLAP_FOUND)
   if (TCLAP_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find TCLAP")
    endif(TCLAP_FIND_REQUIRED)
endif(TCLAP_FOUND) 

set(TCLAP_INCLUDE_DIR "${TCLAP_INCLUDE_DIR}/..")
