# =============================================================================
# add_compare_files_test(
#     NAME        <test_name>
#     TEST_FILE   <path>
#     ORIG_FILE   <path>
#     [ORIG_START <start_line> ORIG_END <end_line>]
#     [MIN_SIZE   <bytes>]
#     [MAX_DIFF   <lines>]
# )
# =============================================================================
function(add_compare_files_test)
  set(options)
  set(oneValueArgs NAME TEST_FILE ORIG_FILE ORIG_START ORIG_END MIN_SIZE MAX_DIFF)
  set(multiValueArgs)
  cmake_parse_arguments(ACFT "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  foreach(VAR NAME TEST_FILE ORIG_FILE)
    if(NOT ACFT_${VAR})
      message(FATAL_ERROR "Argument ${VAR} is mandatory.")
    endif()
  endforeach()

  if((ACFT_ORIG_START AND NOT ACFT_ORIG_END) OR (ACFT_ORIG_END AND NOT ACFT_ORIG_START))
    message(FATAL_ERROR "ORIG_START and ORIG_END must be given together or both omitted.")
  endif()

  set(cmd ${CMAKE_COMMAND}
    -DTEST_FILE=${ACFT_TEST_FILE}
    -DORIG_FILE=${ACFT_ORIG_FILE}
  )

  foreach(VAR ORIG_START ORIG_END MIN_SIZE MAX_DIFF)
    if(ACFT_${VAR})
      list(APPEND cmd -D${VAR}=${ACFT_${VAR}})
    endif()
  endforeach()

  list(APPEND cmd
    -P ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/CompareFiles.cmake
  )

  add_test(NAME ${ACFT_NAME} COMMAND ${cmd})
endfunction()
