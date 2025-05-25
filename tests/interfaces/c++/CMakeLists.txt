
if(NOT DEFINED CEA_CORE)
  message(FATAL_ERROR "This file must be loaded from parent CMakeLists.txt.")
endif()


set(JOB_NAME test-c++-interface)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


set(TEST_EXECUTABLE ${JOB_NAME}_exe)


add_executable(${TEST_EXECUTABLE} EXCLUDE_FROM_ALL)

set(EXAMPLE_DIR "${CMAKE_SOURCE_DIR}/examples/c++")

target_sources(${TEST_EXECUTABLE}
  PRIVATE
  src/main.cpp
  src/legacy_examples.cpp
  ${EXAMPLE_DIR}/src/example-01.cpp
  ${EXAMPLE_DIR}/src/example-02.cpp
  ${EXAMPLE_DIR}/src/example-08.cpp
  ${EXAMPLE_DIR}/src/example-09.cpp
  ${EXAMPLE_DIR}/src/example-10.cpp
  ${EXAMPLE_DIR}/src/example-11.cpp
  ${EXAMPLE_DIR}/src/example-12.cpp
  ${EXAMPLE_DIR}/src/example-13.cpp
)

target_include_directories(${TEST_EXECUTABLE} PRIVATE ${CMAKE_SOURCE_DIR}/src/interfaces/c++ ${EXAMPLE_DIR}/src)

target_link_libraries(${TEST_EXECUTABLE} PRIVATE ${CEA_SHARED_LIBRARY})


add_custom_target(run-${JOB_NAME}
  COMMAND cp -a ${CMAKE_SOURCE_DIR}/tests/orig/cea2.inp .
  COMMAND cp -a ${CMAKE_SOURCE_DIR}/tests/thermp/test-thermp.inp .
  COMMAND cp -a ${CMAKE_SOURCE_DIR}/tests/deton/test-deton.inp .
  COMMAND cp -a ${CMAKE_SOURCE_DIR}/tests/shock/test-shock.inp .
  COMMAND cp -a ${CMAKE_SOURCE_DIR}/tests/rocket/test-rocket.inp .
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
add_test(NAME diff-c++-example-01 COMMAND bash -c "[ $(cat example-01.out | wc -l) -ge 200 ] && [ $(awk 'NR ==    1, NR ==  245' ${ORIG_OUT} | diff -Bb - example-01.out | wc -l) -le 52 ]")
add_test(NAME diff-c++-example-02 COMMAND bash -c "[ $(cat example-02.out | wc -l) -ge 160 ] && [ $(awk 'NR ==  250, NR ==  434' ${ORIG_OUT} | diff -Bb - example-02.out | wc -l) -le 42 ]")
add_test(NAME diff-c++-example-08 COMMAND bash -c "[ $(cat example-08.out | wc -l) -ge 180 ] && [ $(awk 'NR == 1621, NR == 1800' ${ORIG_OUT} | diff -Bb - example-08.out | wc -l) -le 18 ]")
add_test(NAME diff-c++-example-09 COMMAND bash -c "[ $(cat example-09.out | wc -l) -ge 193 ] && [ $(awk 'NR == 1834, NR == 2026' ${ORIG_OUT} | diff -Bb - example-09.out | wc -l) -le 18 ]")
add_test(NAME diff-c++-example-10 COMMAND bash -c "[ $(cat example-10.out | wc -l) -ge 122 ] && [ $(awk 'NR == 2057, NR == 2178' ${ORIG_OUT} | diff -Bb - example-10.out | wc -l) -le 22 ]")
add_test(NAME diff-c++-example-11 COMMAND bash -c "[ $(cat example-11.out | wc -l) -ge 145 ] && [ $(awk 'NR == 2219, NR == 2363' ${ORIG_OUT} | diff -Bb - example-11.out | wc -l) -le 20 ]")
add_test(NAME diff-c++-example-12 COMMAND bash -c "[ $(cat example-12.out | wc -l) -ge 360 ] && [ $(awk 'NR == 2415, NR == 2774' ${ORIG_OUT} | diff -Bb - example-12.out | wc -l) -le 26 ]")
add_test(NAME diff-c++-example-13 COMMAND bash -c "[ $(cat example-13.out | wc -l) -ge 185 ] && [ $(awk 'NR == 2819, NR == 3003' ${ORIG_OUT} | diff -Bb - example-13.out | wc -l) -le 22 ]")

foreach(i IN ITEMS 01 02 08 09 10 11 12 13)
  set_tests_properties(diff-c++-example-${i} PROPERTIES FIXTURES_REQUIRED ${JOB_NAME}_fixtures)
endforeach()

add_test(
  NAME diff-c++-example-12-plt
  COMMAND diff -Bb ${CMAKE_SOURCE_DIR}/tests/orig/cea2.plt example-12.plt
)
set_tests_properties(diff-c++-example-12-plt PROPERTIES FIXTURES_REQUIRED ${JOB_NAME}_fixtures)


foreach(TEST_PATH IN ITEMS "orig/cea2" "thermp/test-thermp" "deton/test-deton" "shock/test-shock" "rocket/test-rocket")
  get_filename_component(JOB_NAME ${TEST_PATH} NAME)
  foreach(EXT IN ITEMS "out" "plt")
    if(EXISTS ${CMAKE_SOURCE_DIR}/tests/${TEST_PATH}.${EXT})
      add_test(
        NAME c++/${JOB_NAME}.${EXT}
        COMMAND diff -Bb ${CMAKE_SOURCE_DIR}/tests/${TEST_PATH}.${EXT} ${CMAKE_CURRENT_BINARY_DIR}/${JOB_NAME}.${EXT}
      )
      set_tests_properties(c++/${JOB_NAME}.${EXT} PROPERTIES FIXTURES_REQUIRED ${JOB_NAME}_fixtures)
    endif()
  endforeach()
endforeach()
