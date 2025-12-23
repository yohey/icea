# =============================================================================
# cmake -P CompareFiles.cmake
#       -DTEST_FILE=<path>
#       -DORIG_FILE=<path>
#       [-DORIG_START=<start_line> -DORIG_END=<end_line>]
#       [-DMIN_SIZE=<bytes>]
#       [-DMAX_DIFF=<lines>]
# =============================================================================
cmake_minimum_required(VERSION 2.6)

# --- Required arguments check ------------------------------------------------
if(NOT DEFINED TEST_FILE OR NOT DEFINED ORIG_FILE)
  message(FATAL_ERROR "TEST_FILE and ORIG_FILE must be defined.")
endif()
if(NOT EXISTS "${TEST_FILE}")
  message(FATAL_ERROR "TEST_FILE '${TEST_FILE}' does not exist.")
endif()
if(NOT EXISTS "${ORIG_FILE}")
  message(FATAL_ERROR "ORIG_FILE '${ORIG_FILE}' does not exist.")
endif()

set(_ORIG_FILE_NAME "${ORIG_FILE}")

# --- MIN_SIZE check ----------------------------------------------------------
if(DEFINED MIN_SIZE)
  file(SIZE "${TEST_FILE}" _size)
  if(_size LESS MIN_SIZE)
    message(FATAL_ERROR "File size ${_size} B is smaller than MIN_SIZE = ${MIN_SIZE} B.")
  endif()
endif()

# --- ORIG_START / ORIG_END extraction ------------------------------------------
if(DEFINED ORIG_START AND DEFINED ORIG_END)
  math(EXPR _slice_len "${ORIG_END} - ${ORIG_START} + 1")
  if(_slice_len LESS_EQUAL 0)
    message(FATAL_ERROR "ORIG_END must be ≥ ORIG_START.")
  endif()

  file(STRINGS "${ORIG_FILE}" _orig_lines)
  math(EXPR _start_idx "${ORIG_START} - 1")
  list(LENGTH _orig_lines _orig_len)
  if(ORIG_END GREATER _orig_len)
    message(FATAL_ERROR "ORIG_END (${ORIG_END}) exceeds total lines in ORIG_FILE (${_orig_len}).")
  endif()
  list(SUBLIST _orig_lines ${_start_idx} ${_slice_len} _orig_lines)

  set(_orig_slice "${CMAKE_BINARY_DIR}/tmp_orig_slice.txt")
  file(WRITE "${_orig_slice}" "")
  foreach(line IN LISTS _orig_lines)
    file(APPEND "${_orig_slice}" "${line}\n")
  endforeach()
  set(ORIG_FILE "${_orig_slice}")
elseif(DEFINED ORIG_START OR DEFINED ORIG_END)
  message(FATAL_ERROR "ORIG_START and ORIG_END must be given together.")
endif()

# --- Perform comparison ------------------------------------------------------
if(WIN32)
  execute_process(
    COMMAND powershell -NoProfile -ExecutionPolicy Bypass
                       -File "${CMAKE_CURRENT_LIST_DIR}/Compare-Files.ps1"
                       -TestFile "${TEST_FILE}" -OrigFile "${ORIG_FILE}"
    RESULT_VARIABLE _diff_ret
    OUTPUT_VARIABLE _diff_out
    ERROR_VARIABLE  _diff_err
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE)
else()
  execute_process(
    COMMAND diff -b -B "${ORIG_FILE}" "${TEST_FILE}"
    RESULT_VARIABLE _diff_ret
    OUTPUT_VARIABLE _diff_out
    ERROR_VARIABLE  _diff_err
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE)
endif()

# --- MAX_DIFF judgment or full match check -----------------------------------
if(DEFINED MAX_DIFF)
  string(REGEX MATCHALL "\n" _new_lines "${_diff_out}\n")
  list(LENGTH _new_lines _diff_cnt)
  message(STATUS "_diff_cnt = ${_diff_cnt}")

  if(_diff_cnt GREATER MAX_DIFF)
    message(STATUS "Diff lines (${_diff_cnt}) exceed MAX_DIFF = ${MAX_DIFF}:")
    message(STATUS "${_diff_out}")
    if(_diff_err)
      message(STATUS "${_diff_err}")
    endif()
    message(FATAL_ERROR "File difference exceeds MAX_DIFF.")
  endif()
elseif(NOT _diff_ret EQUAL 0)
  message(STATUS "Files differ:\n${_diff_out}")
  if(_diff_err)
    message(STATUS "${_diff_err}")
  endif()
  if(DEFINED ORIG_START)
    message(FATAL_ERROR "Files differ: '${TEST_FILE}' ≠ '${_ORIG_FILE_NAME}' (ll. ${ORIG_START}--${ORIG_END}).")
  else()
    message(FATAL_ERROR "Files differ: '${TEST_FILE}' ≠ '${_ORIG_FILE_NAME}'.")
  endif()
endif()
