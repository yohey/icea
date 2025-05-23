
set(JOB_NAME test-rocket)

add_custom_target(run-${JOB_NAME}
  COMMAND cp -a ${CMAKE_CURRENT_SOURCE_DIR}/${JOB_NAME}.inp .
  COMMAND echo ${JOB_NAME} | ${CMAKE_BINARY_DIR}/${CEA_EXECUTABLE} --legacy-mode --thermo-lib ${CMAKE_BINARY_DIR}/thermo.lib --trans-lib ${CMAKE_BINARY_DIR}/trans.lib
  DEPENDS ${CEA_EXECUTABLE}
)

add_test(
  NAME ${JOB_NAME}
  COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config "$<CONFIG>" --target run-${JOB_NAME}
)

set_tests_properties(${JOB_NAME} PROPERTIES FIXTURES_SETUP ${JOB_NAME}_fixtures)

foreach(EXT IN ITEMS "out" "plt")
  add_test(
    NAME ${JOB_NAME}.${EXT}
    COMMAND diff -Bb ${CMAKE_CURRENT_SOURCE_DIR}/${JOB_NAME}.${EXT} ${CMAKE_CURRENT_BINARY_DIR}/${JOB_NAME}.${EXT}
  )
  set_tests_properties(${JOB_NAME}.${EXT} PROPERTIES FIXTURES_REQUIRED ${JOB_NAME}_fixtures)
endforeach()
