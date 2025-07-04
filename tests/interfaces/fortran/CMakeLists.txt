
if(NOT DEFINED CEA_CORE)
  message(FATAL_ERROR "This file must be loaded from parent CMakeLists.txt.")
endif()

include(${CMAKE_SOURCE_DIR}/tests/cmake/AddCompareFilesTest.cmake)


set(JOB_NAME test-fortran-interface)

set(TEST_EXECUTABLE ${JOB_NAME}_exe)

add_executable(${TEST_EXECUTABLE} EXCLUDE_FROM_ALL)

set(EXAMPLE_DIR "${CMAKE_SOURCE_DIR}/examples/fortran")

target_sources(${TEST_EXECUTABLE}
  PRIVATE
  ${EXAMPLE_DIR}/src/main.f90
  ${EXAMPLE_DIR}/src/example-01.f90
  ${EXAMPLE_DIR}/src/example-02.f90
  ${EXAMPLE_DIR}/src/example-03.f90
  ${EXAMPLE_DIR}/src/example-04.f90
  ${EXAMPLE_DIR}/src/example-05.f90
  ${EXAMPLE_DIR}/src/example-06.f90
  ${EXAMPLE_DIR}/src/example-07.f90
  ${EXAMPLE_DIR}/src/example-08.f90
  ${EXAMPLE_DIR}/src/example-09.f90
  ${EXAMPLE_DIR}/src/example-10.f90
  ${EXAMPLE_DIR}/src/example-11.f90
  ${EXAMPLE_DIR}/src/example-12.f90
  ${EXAMPLE_DIR}/src/example-13.f90
  ${EXAMPLE_DIR}/src/example-14.f90
)

target_compile_options(${TEST_EXECUTABLE}
  PRIVATE
  ${CEA_Fortran_FLAGS}
)

target_include_directories(${TEST_EXECUTABLE} PRIVATE ${CMAKE_BINARY_DIR})

target_link_libraries(${TEST_EXECUTABLE} PRIVATE ${CEA_CORE})


add_custom_target(run-${JOB_NAME}
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/thermo.lib ./thermo.lib
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/trans.lib  ./trans.lib
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${TEST_EXECUTABLE}
  DEPENDS ${TEST_EXECUTABLE}
)


add_test(
  NAME ${JOB_NAME}
  COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config "$<CONFIG>" --target run-${JOB_NAME}
)

set_tests_properties(${JOB_NAME} PROPERTIES FIXTURES_SETUP ${JOB_NAME}_fixtures)


set(ORIG_OUT "${CMAKE_SOURCE_DIR}/tests/orig/cea2.out")
add_compare_files_test(NAME diff-fortran-example-01 TEST_FILE example-01.out ORIG_FILE ${ORIG_OUT} ORIG_START   53 ORIG_END  245 MIN_SIZE  9318 MAX_DIFF 18)
add_compare_files_test(NAME diff-fortran-example-02 TEST_FILE example-02.out ORIG_FILE ${ORIG_OUT} ORIG_START  279 ORIG_END  434 MIN_SIZE  5692 MAX_DIFF 14)
add_compare_files_test(NAME diff-fortran-example-03 TEST_FILE example-03.out ORIG_FILE ${ORIG_OUT} ORIG_START  491 ORIG_END  654 MIN_SIZE  8168 MAX_DIFF 14)
add_compare_files_test(NAME diff-fortran-example-04 TEST_FILE example-04.out ORIG_FILE ${ORIG_OUT} ORIG_START  711 ORIG_END  870 MIN_SIZE  7022 MAX_DIFF 14)
add_compare_files_test(NAME diff-fortran-example-05 TEST_FILE example-05.out ORIG_FILE ${ORIG_OUT} ORIG_START  937 ORIG_END 1183 MIN_SIZE 15074 MAX_DIFF 14)
add_compare_files_test(NAME diff-fortran-example-06 TEST_FILE example-06.out ORIG_FILE ${ORIG_OUT} ORIG_START 1214 ORIG_END 1363 MIN_SIZE  5835 MAX_DIFF 14)
add_compare_files_test(NAME diff-fortran-example-07 TEST_FILE example-07.out ORIG_FILE ${ORIG_OUT} ORIG_START 1402 ORIG_END 1581 MIN_SIZE  9615 MAX_DIFF 18)
add_compare_files_test(NAME diff-fortran-example-08 TEST_FILE example-08.out ORIG_FILE ${ORIG_OUT} ORIG_START 1621 ORIG_END 1800 MIN_SIZE  8004 MAX_DIFF 18)
add_compare_files_test(NAME diff-fortran-example-09 TEST_FILE example-09.out ORIG_FILE ${ORIG_OUT} ORIG_START 1834 ORIG_END 2026 MIN_SIZE  8788 MAX_DIFF 18)
add_compare_files_test(NAME diff-fortran-example-10 TEST_FILE example-10.out ORIG_FILE ${ORIG_OUT} ORIG_START 2057 ORIG_END 2178 MIN_SIZE  6100 MAX_DIFF 22)
add_compare_files_test(NAME diff-fortran-example-11 TEST_FILE example-11.out ORIG_FILE ${ORIG_OUT} ORIG_START 2219 ORIG_END 2363 MIN_SIZE  7177 MAX_DIFF 20)
add_compare_files_test(NAME diff-fortran-example-12 TEST_FILE example-12.out ORIG_FILE ${ORIG_OUT} ORIG_START 2415 ORIG_END 2774 MIN_SIZE 17022 MAX_DIFF 26)
add_compare_files_test(NAME diff-fortran-example-13 TEST_FILE example-13.out ORIG_FILE ${ORIG_OUT} ORIG_START 2819 ORIG_END 3003 MIN_SIZE  9590 MAX_DIFF 22)
add_compare_files_test(NAME diff-fortran-example-14 TEST_FILE example-14.out ORIG_FILE ${ORIG_OUT} ORIG_START 3043 ORIG_END 3367 MIN_SIZE 19036 MAX_DIFF 14)

foreach(i IN ITEMS 01 02 03 04 05 06 07 08 09 10 11 12 13 14)
  set_tests_properties(diff-fortran-example-${i} PROPERTIES FIXTURES_REQUIRED ${JOB_NAME}_fixtures)
endforeach()

add_compare_files_test(
  NAME      diff-fortran-example-12-plt
  TEST_FILE example-12.plt
  ORIG_FILE ${CMAKE_SOURCE_DIR}/tests/orig/cea2.plt
)
set_tests_properties(diff-fortran-example-12-plt PROPERTIES FIXTURES_REQUIRED ${JOB_NAME}_fixtures)
