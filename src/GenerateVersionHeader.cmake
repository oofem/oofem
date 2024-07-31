# Credits to https://github.com/nocnokneo/cmake-git-versioning-example
# Copyright (c) 2021 Taylor Braun-Jones
if(GIT_EXECUTABLE)
  #get_filename_component(SRC_DIR ${CMAKE_SOURCE_DIR} DIRECTORY)
  message (${OOFEM_SRC_DIR})
  execute_process(
    #COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
    COMMAND ${GIT_EXECUTABLE} describe --tags --dirty
    WORKING_DIRECTORY ${OOFEM_SRC_DIR}
    OUTPUT_VARIABLE GIT_HASH
    RESULT_VARIABLE GIT_HASH_ERROR_CODE
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT GIT_HASH_ERROR_CODE)
    set(OOFEM_GIT_HASH ${GIT_HASH})
  endif()

  # get repository URL & branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} remote get-url origin
    WORKING_DIRECTORY ${OOFEM_SRC_DIR}
    OUTPUT_VARIABLE GIT_REPOURL
    RESULT_VARIABLE GIT_REPOURL_ERROR_CODE
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT GIT_REPOURL_ERROR_CODE)
    set(OOFEM_GIT_REPOURL ${GIT_REPOURL})
  endif()
  
  execute_process(
    COMMAND ${GIT_EXECUTABLE} branch --show-current
    WORKING_DIRECTORY ${OOFEM_SRC_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    RESULT_VARIABLE GIT_BRANCHL_ERROR_CODE
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT GIT_BRANCH_ERROR_CODE)
    set(OOFEM_GIT_BRANCH ${GIT_BRANCH})
  endif()



endif()
# Final fallback: Just use a bogus version string that is semantically older
# than anything else and spit out a warning to the developer.

if(NOT DEFINED OOFEM_GIT_HASH)
  set(OOFEM_GIT_HASH unknown)
  message(WARNING "Failed to determine OOFEM_GIT_HASH from Git log. Using default hash \"${OOFEM_GIT_HASH}\".")
endif()

if(NOT DEFINED OOFEM_GIT_REPOURL)
  set(OOFEM_GIT_REPOURL unknown)
  message(WARNING "Failed to determine OOFEM_GIT_REPOURL from Git log.")
endif()

if(NOT DEFINED OOFEM_GIT_BRANCH)
  set(OOFEM_GIT_BRANCH unknown)
  message(WARNING "Failed to determine OOFEM_GIT_HASH from Git log.")
endif()
configure_file(${SRC} ${DST} @ONLY)