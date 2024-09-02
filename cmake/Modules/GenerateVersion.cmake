#######################################################################
# generate git-version signature
#######################################################################
include_directories(${CMAKE_BINARY_DIR})
find_package(Git)
add_custom_target(version
  ${CMAKE_COMMAND} -D SRC=${CMAKE_SOURCE_DIR}/src/oofem_version.h.in
                   -D DST=${CMAKE_BINARY_DIR}/oofem_version.h
                   -D OOFEM_SRC_DIR=${CMAKE_SOURCE_DIR}
                   -D GIT_EXECUTABLE=${GIT_EXECUTABLE}
                   -P ${CMAKE_SOURCE_DIR}/src/GenerateVersionHeader.cmake
  )