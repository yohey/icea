
set(JOB_NAME test-c++)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


set(TEST_EXECUTABLE ${JOB_NAME})


add_executable(${TEST_EXECUTABLE} EXCLUDE_FROM_ALL)

target_sources(${TEST_EXECUTABLE}
  PRIVATE
  src/main.cpp
)

target_include_directories(${TEST_EXECUTABLE} PRIVATE ${CMAKE_SOURCE_DIR}/src/interfaces/c++)

target_link_libraries(${TEST_EXECUTABLE} PRIVATE ${CEA_CXX_LIBRARY})


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
