cmake_minimum_required(VERSION 3.5)

enable_language(Fortran)
project(CEA Fortran)


set(CEA_EXECUTABLE cea2)


add_executable(${CEA_EXECUTABLE})

target_sources(${CEA_EXECUTABLE}
  PRIVATE
  src/cea2.f90
  src/mod_cea.f90
  src/mod_general.f90
  src/mod_legacy_io.f90
)

target_compile_options(${CEA_EXECUTABLE}
  PRIVATE
  "-std=f2023"
)


foreach(LIB IN ITEMS "thermo" "trans")
  add_custom_target(${LIB} ALL
    COMMAND echo ${LIB} | ./${CEA_EXECUTABLE} > /dev/null
    DEPENDS ${CEA_EXECUTABLE} ${LIB}.inp
  )

  add_custom_command(
    OUTPUT ${LIB}.inp
    COMMAND cp -a ${CMAKE_SOURCE_DIR}/src/${LIB}.inp ${CMAKE_BINARY_DIR}/
  )
endforeach()


enable_testing()

foreach(INP IN ITEMS "thermo" "trans")
  foreach(EXT IN ITEMS "lib" "out")
    add_test(
      NAME ${INP}.${EXT}
      COMMAND diff ${CMAKE_SOURCE_DIR}/test/${INP}.${EXT} ${CMAKE_BINARY_DIR}/${INP}.${EXT}
    )
  endforeach()
endforeach()

add_custom_target(run-cea2
  COMMAND echo cea2 | ./${CEA_EXECUTABLE} > /dev/null
  DEPENDS ${CEA_EXECUTABLE} cea2.inp
)

add_custom_command(
  OUTPUT cea2.inp
  COMMAND cp -a ${CMAKE_SOURCE_DIR}/test/cea2.inp ${CMAKE_BINARY_DIR}/
)

add_test(
  NAME cea2
  COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config "$<CONFIG>" --target run-cea2
)

set_tests_properties(cea2 PROPERTIES FIXTURES_SETUP test_fixtures)

foreach(EXT IN ITEMS "out" "plt")
  add_test(
    NAME cea2.${EXT}
    COMMAND diff -B ${CMAKE_SOURCE_DIR}/test/cea2.${EXT} ${CMAKE_BINARY_DIR}/cea2.${EXT}
  )
  set_tests_properties(cea2.${EXT} PROPERTIES FIXTURES_REQUIRED test_fixtures)
endforeach()


install(TARGETS ${CEA_EXECUTABLE} RUNTIME DESTINATION bin)
