
if(NOT DEFINED CEA_CORE)
  message(FATAL_ERROR "This file must be loaded from parent CMakeLists.txt.")
endif()

include(${CMAKE_SOURCE_DIR}/tests/cmake/AddCompareFilesTest.cmake)

foreach(INP IN ITEMS "thermo" "trans")
  foreach(EXT IN ITEMS "lib" "out")
    add_compare_files_test(
      NAME      ${INP}.${EXT}
      TEST_FILE ${CMAKE_BINARY_DIR}/${INP}.${EXT}
      ORIG_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${INP}.${EXT}
    )
  endforeach()
endforeach()

set(JOB_NAME cea2)

add_custom_target(run-${JOB_NAME}
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${JOB_NAME}.inp .
  COMMAND echo ${JOB_NAME} | \"${CMAKE_BINARY_DIR}/${CEA_EXECUTABLE}\" --legacy-mode --thermo-lib ${CMAKE_BINARY_DIR}/thermo.lib --trans-lib ${CMAKE_BINARY_DIR}/trans.lib
  DEPENDS ${CEA_EXECUTABLE}
)

add_test(
  NAME ${JOB_NAME}
  COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config "$<CONFIG>" --target run-${JOB_NAME}
)

set_tests_properties(${JOB_NAME} PROPERTIES FIXTURES_SETUP ${JOB_NAME}_fixtures)

foreach(EXT IN ITEMS "out" "plt")
  add_compare_files_test(
    NAME      ${JOB_NAME}.${EXT}
    TEST_FILE ${CMAKE_CURRENT_BINARY_DIR}/${JOB_NAME}.${EXT}
    ORIG_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${JOB_NAME}.${EXT}
  )
  set_tests_properties(${JOB_NAME}.${EXT} PROPERTIES FIXTURES_REQUIRED ${JOB_NAME}_fixtures)
endforeach()
