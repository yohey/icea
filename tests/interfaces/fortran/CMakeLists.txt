
if(NOT DEFINED CEA_CORE)
  message(FATAL_ERROR "This file must be loaded from parent CMakeLists.txt.")
endif()


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
  COMMAND ln -fs ${CMAKE_BINARY_DIR}/thermo.lib .
  COMMAND ln -fs ${CMAKE_BINARY_DIR}/trans.lib .
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${TEST_EXECUTABLE}
  DEPENDS ${TEST_EXECUTABLE}
)


add_test(
  NAME ${JOB_NAME}
  COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config "$<CONFIG>" --target run-${JOB_NAME}
)

set_tests_properties(${JOB_NAME} PROPERTIES FIXTURES_SETUP ${JOB_NAME}_fixtures)


set(ORIG_OUT "${CMAKE_SOURCE_DIR}/tests/orig/cea2.out")
add_test(NAME diff-fortran-example-01 COMMAND bash -c "[ $(cat example-01.out | wc -l) -ge 201 ] && [ $(awk 'NR ==   53, NR ==  245' ${ORIG_OUT} | diff -Bb - example-01.out | wc -l) -le 18 ]")
add_test(NAME diff-fortran-example-02 COMMAND bash -c "[ $(cat example-02.out | wc -l) -ge 164 ] && [ $(awk 'NR ==  279, NR ==  434' ${ORIG_OUT} | diff -Bb - example-02.out | wc -l) -le 14 ]")
add_test(NAME diff-fortran-example-03 COMMAND bash -c "[ $(cat example-03.out | wc -l) -ge 172 ] && [ $(awk 'NR ==  491, NR ==  654' ${ORIG_OUT} | diff -Bb - example-03.out | wc -l) -le 14 ]")
add_test(NAME diff-fortran-example-04 COMMAND bash -c "[ $(cat example-04.out | wc -l) -ge 168 ] && [ $(awk 'NR ==  711, NR ==  870' ${ORIG_OUT} | diff -Bb - example-04.out | wc -l) -le 14 ]")
add_test(NAME diff-fortran-example-05 COMMAND bash -c "[ $(cat example-05.out | wc -l) -ge 255 ] && [ $(awk 'NR ==  937, NR == 1183' ${ORIG_OUT} | diff -Bb - example-05.out | wc -l) -le 14 ]")
add_test(NAME diff-fortran-example-06 COMMAND bash -c "[ $(cat example-06.out | wc -l) -ge 158 ] && [ $(awk 'NR == 1214, NR == 1363' ${ORIG_OUT} | diff -Bb - example-06.out | wc -l) -le 14 ]")
add_test(NAME diff-fortran-example-07 COMMAND bash -c "[ $(cat example-07.out | wc -l) -ge 188 ] && [ $(awk 'NR == 1402, NR == 1581' ${ORIG_OUT} | diff -Bb - example-07.out | wc -l) -le 18 ]")
add_test(NAME diff-fortran-example-08 COMMAND bash -c "[ $(cat example-08.out | wc -l) -ge 189 ] && [ $(awk 'NR == 1621, NR == 1800' ${ORIG_OUT} | diff -Bb - example-08.out | wc -l) -le 18 ]")
add_test(NAME diff-fortran-example-09 COMMAND bash -c "[ $(cat example-09.out | wc -l) -ge 202 ] && [ $(awk 'NR == 1834, NR == 2026' ${ORIG_OUT} | diff -Bb - example-09.out | wc -l) -le 18 ]")
add_test(NAME diff-fortran-example-10 COMMAND bash -c "[ $(cat example-10.out | wc -l) -ge 135 ] && [ $(awk 'NR == 2057, NR == 2178' ${ORIG_OUT} | diff -Bb - example-10.out | wc -l) -le 22 ]")
add_test(NAME diff-fortran-example-11 COMMAND bash -c "[ $(cat example-11.out | wc -l) -ge 154 ] && [ $(awk 'NR == 2219, NR == 2363' ${ORIG_OUT} | diff -Bb - example-11.out | wc -l) -le 20 ]")
add_test(NAME diff-fortran-example-12 COMMAND bash -c "[ $(cat example-12.out | wc -l) -ge 369 ] && [ $(awk 'NR == 2415, NR == 2774' ${ORIG_OUT} | diff -Bb - example-12.out | wc -l) -le 26 ]")
add_test(NAME diff-fortran-example-13 COMMAND bash -c "[ $(cat example-13.out | wc -l) -ge 194 ] && [ $(awk 'NR == 2819, NR == 3003' ${ORIG_OUT} | diff -Bb - example-13.out | wc -l) -le 22 ]")
add_test(NAME diff-fortran-example-14 COMMAND bash -c "[ $(cat example-14.out | wc -l) -ge 333 ] && [ $(awk 'NR == 3043, NR == 3367' ${ORIG_OUT} | diff -Bb - example-14.out | wc -l) -le 14 ]")

foreach(i IN ITEMS 01 02 03 04 05 06 07 08 09 10 11 12 13 14)
  set_tests_properties(diff-fortran-example-${i} PROPERTIES FIXTURES_REQUIRED ${JOB_NAME}_fixtures)
endforeach()

add_test(
  NAME diff-fortran-example-12-plt
  COMMAND diff -Bb ${CMAKE_SOURCE_DIR}/tests/orig/cea2.plt example-12.plt
)
set_tests_properties(diff-fortran-example-12-plt PROPERTIES FIXTURES_REQUIRED ${JOB_NAME}_fixtures)
