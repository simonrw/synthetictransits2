#To use this library you must add it as a subdirectory, then use 
#
#  MODELGEN_LIBRARIES
#  MODELGEN_INCLUDE_DIR
# 
# to get the library and include files respectively linked up. 
# 
# Also you must add the MODELTARGET name to the dependency list
# so it compiles before any code that uses it.


# Name of the library to create
set(TARGET Modelgen)

##############################################################
# YOU SHOULD NOT NEED TO CHANGE ANYTHING BELOW HERE
##############################################################

# Need to create this variable in the parent scope too
set(MODELTARGET ${TARGET} PARENT_SCOPE)
project(${TARGET})

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    )

add_library(${TARGET}
    STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/GenerateModel.cpp
    )

# set up some nice variables for the parent scope
set(MODELGEN_LIBRARIES ${${TARGET}_BINARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}${TARGET}${CMAKE_STATIC_LIBRARY_SUFFIX} PARENT_SCOPE)
set(MODELGEN_INCLUDE_DIR ${${TARGET}_SOURCE_DIR}/include PARENT_SCOPE)
