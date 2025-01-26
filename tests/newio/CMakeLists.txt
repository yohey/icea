
set(JOB_NAME test-newio)

set(TEST_EXECUTABLE ${JOB_NAME}_exe)


add_executable(${TEST_EXECUTABLE} EXCLUDE_FROM_ALL)

target_sources(${TEST_EXECUTABLE}
  PRIVATE
  test_newio.f90
)

target_compile_options(${TEST_EXECUTABLE}
  PRIVATE
  "-cpp"
  "-std=f2023"
)

target_include_directories(${TEST_EXECUTABLE} PRIVATE ${CMAKE_BINARY_DIR})

target_link_libraries(${TEST_EXECUTABLE} PRIVATE ${CEA_LIBRARY})


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
