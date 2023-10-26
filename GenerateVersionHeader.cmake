# Credits to https://github.com/nocnokneo/cmake-git-versioning-example
# Copyright (c) 2021 Taylor Braun-Jones
if(GIT_EXECUTABLE)
  get_filename_component(SRC_DIR ${SRC} DIRECTORY)
  # Generate a git-describe version string from Git repository tags
  #execute_process(
  #  COMMAND ${GIT_EXECUTABLE} describe --tags --dirty --match "v*"
  #  WORKING_DIRECTORY ${SRC_DIR}
  #  OUTPUT_VARIABLE GIT_DESCRIBE_VERSION
  #  RESULT_VARIABLE GIT_DESCRIBE_ERROR_CODE
  #  OUTPUT_STRIP_TRAILING_WHITESPACE
  #  )
  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
    WORKING_DIRECTORY ${SRC_DIR}
    OUTPUT_VARIABLE GIT_HASH
    RESULT_VARIABLE GIT_HASH_ERROR_CODE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  if(NOT GIT_HASH_ERROR_CODE)
    set(OOFEM_GIT_HASH ${GIT_HASH})
  endif()
endif()

# Final fallback: Just use a bogus version string that is semantically older
# than anything else and spit out a warning to the developer.
if(NOT DEFINED OOFEM_GIT_HASH)
  set(OOFEM_GIT_HASH unknown)
  message(WARNING "Failed to determine OOFEM_GIT_HASH from Git log. Using default hash \"${OOFEM_GIT_HASH}\".")
endif()

configure_file(${SRC} ${DST} @ONLY)