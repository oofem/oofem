
find_package(Git)
if(SKBUILD_PROJECT_VERSION)
  message(STATUS "Using SKBUILD_PROJECT_VERSION = ${SKBUILD_PROJECT_VERSION} for version information")
  string(REGEX MATCHALL "^([0-9]*)\.([0-9]*)(\.([0-9]*)\.?(.*))?$" _matches ${SKBUILD_PROJECT_VERSION})
  set(OOFEM_VERSION_MAJOR ${CMAKE_MATCH_1})
  set(OOFEM_VERSION_MINOR ${CMAKE_MATCH_2})
  set(OOFEM_VERSION_PATCH ${CMAKE_MATCH_4})
  set(OOFEM_GIT_HASH "skbuild")
  set(OOFEM_GIT_REPOURL "python-sdist")
  set(OOFEM_GIT_BRANCH "python-sdist")
elseif(GIT_EXECUTABLE)
  message(STATUS "Using git describe for version information")
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --dirty
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    RESULT_VARIABLE GIT_HASH_ERROR_CODE
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT GIT_HASH_ERROR_CODE)
    set(OOFEM_GIT_HASH ${GIT_HASH})
    string(REGEX MATCHALL "^v([0-9]*)\.([0-9]*)\.*([0-9]*)" _matches ${GIT_HASH})
    set(OOFEM_VERSION_MAJOR ${CMAKE_MATCH_1})
    set(OOFEM_VERSION_MINOR ${CMAKE_MATCH_2})
    if (NOT CMAKE_MATCH_3)
      set(OOFEM_VERSION_PATCH 0)
    else()
      set(OOFEM_VERSION_PATCH ${CMAKE_MATCH_3})
    endif()
  endif()

  # get repository URL & branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} remote get-url origin
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_REPOURL
    RESULT_VARIABLE GIT_REPOURL_ERROR_CODE
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT GIT_REPOURL_ERROR_CODE)
    set(OOFEM_GIT_REPOURL ${GIT_REPOURL})
  endif()
  
  execute_process(
    COMMAND ${GIT_EXECUTABLE} branch --show-current
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    RESULT_VARIABLE GIT_BRANCH_ERROR_CODE
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT GIT_BRANCH_ERROR_CODE)
    set(OOFEM_GIT_BRANCH ${GIT_BRANCH})
  endif()
endif()
if(NOT DEFINED OOFEM_GIT_HASH)
  message(WARNING "Using <unknown> for version information (git failed or not found)")
  set(OOFEM_VERSION_MAJOR "<?>")
  set(OOFEM_VERSION_MINOR "<?>")
  set(OOFEM_VERSION_PATCH "<?>")
  set(OOFEM_GIT_HASH "<unknown>")
  set(OOFEM_GIT_REPOURL "<unknown>")
  set(OOFEM_GIT_BRANCH "<unknown>")
endif()
message(STATUS "Versions are: ${OOFEM_GIT_HASH} (${OOFEM_VERSION_MAJOR}.${OOFEM_VERSION_MINOR}.${OOFEM_VERSION_PATCH}), repo: ${OOFEM_GIT_REPOURL}, branch: ${OOFEM_GIT_BRANCH}")

