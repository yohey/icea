
if(NOT DEFINED CEA_CORE)
  message(FATAL_ERROR "This file must be loaded from parent CMakeLists.txt.")
endif()


set(JOB_NAME test-python-interface)

set(EXAMPLE_DIR "${CMAKE_SOURCE_DIR}/examples/python")

set("${JOB_NAME}_ENVIRON_LIST"
  "DYLD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/src/interfaces/c++"
  "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
)

include(FindPython)


add_custom_target(copy-${JOB_NAME}-files
  COMMAND ln -s ${CMAKE_BINARY_DIR}/thermo.lib ${CMAKE_CURRENT_BINARY_DIR}/
  COMMAND ln -s ${CMAKE_BINARY_DIR}/trans.lib  ${CMAKE_CURRENT_BINARY_DIR}/
  COMMAND cp -a ${CMAKE_SOURCE_DIR}/src/interfaces/python/icea ${CMAKE_CURRENT_BINARY_DIR}/
)

add_test(
  NAME prepare-${JOB_NAME}
  COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config "$<CONFIG>" --target copy-${JOB_NAME}-files
)


foreach(i IN ITEMS 08 09 10 11 12 13)
  add_test(NAME test-python-example-${i} COMMAND ${Python_INTERPRETER} ${EXAMPLE_DIR}/example-${i}.py)
  set_tests_properties(test-python-example-${i} PROPERTIES ENVIRONMENT "${${JOB_NAME}_ENVIRON_LIST}" DEPENDS prepare-${JOB_NAME})
endforeach()


set(ORIG_OUT "${CMAKE_SOURCE_DIR}/tests/orig/cea2.out")
add_test(NAME diff-python-example-08 COMMAND bash -c "[ $(cat example-08.out | wc -l) -ge 180 ] && [ $(awk 'NR == 1621, NR == 1800' ${ORIG_OUT} | diff -Bb - example-08.out | wc -l) -le 18 ]")
add_test(NAME diff-python-example-09 COMMAND bash -c "[ $(cat example-09.out | wc -l) -ge 193 ] && [ $(awk 'NR == 1834, NR == 2026' ${ORIG_OUT} | diff -Bb - example-09.out | wc -l) -le 18 ]")
add_test(NAME diff-python-example-10 COMMAND bash -c "[ $(cat example-10.out | wc -l) -ge 122 ] && [ $(awk 'NR == 2057, NR == 2178' ${ORIG_OUT} | diff -Bb - example-10.out | wc -l) -le 22 ]")
add_test(NAME diff-python-example-11 COMMAND bash -c "[ $(cat example-11.out | wc -l) -ge 145 ] && [ $(awk 'NR == 2219, NR == 2363' ${ORIG_OUT} | diff -Bb - example-11.out | wc -l) -le 20 ]")
add_test(NAME diff-python-example-12 COMMAND bash -c "[ $(cat example-12.out | wc -l) -ge 360 ] && [ $(awk 'NR == 2415, NR == 2774' ${ORIG_OUT} | diff -Bb - example-12.out | wc -l) -le 26 ]")
add_test(NAME diff-python-example-13 COMMAND bash -c "[ $(cat example-13.out | wc -l) -ge 185 ] && [ $(awk 'NR == 2819, NR == 3003' ${ORIG_OUT} | diff -Bb - example-13.out | wc -l) -le 22 ]")

foreach(i IN ITEMS 08 09 10 11 12 13)
  set_tests_properties(diff-python-example-${i} PROPERTIES DEPENDS test-python-example-${i})
endforeach()

add_test(
  NAME diff-python-example-12-plt
  COMMAND diff -Bb ${CMAKE_SOURCE_DIR}/tests/orig/cea2.plt example-12.plt
)
set_tests_properties(diff-python-example-12-plt PROPERTIES DEPENDS test-python-example-12)
