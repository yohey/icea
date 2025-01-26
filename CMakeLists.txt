cmake_minimum_required(VERSION 3.5)

enable_language(Fortran)
project(CEA Fortran)


set(CEA_LIBRARY libcea)
set(CEA_EXECUTABLE cea)


add_library(${CEA_LIBRARY})

set_target_properties(${CEA_LIBRARY}
  PROPERTIES
  OUTPUT_NAME cea
)

target_sources(${CEA_LIBRARY}
  PRIVATE
  src/mod_cea.f90
  src/mod_types.f90
  src/mod_constants.f90
  src/mod_general.f90
  src/mod_functions.f90
  src/mod_io.f90
  src/mod_legacy_io.f90
)

target_compile_options(${CEA_LIBRARY}
  PRIVATE
  "-cpp"
  "-std=f2023"
)


add_executable(${CEA_EXECUTABLE})

target_sources(${CEA_EXECUTABLE}
  PRIVATE
  src/main.f90
)

target_compile_options(${CEA_EXECUTABLE}
  PRIVATE
  "-cpp"
  "-std=f2023"
)

target_link_libraries(${CEA_EXECUTABLE} PRIVATE ${CEA_LIBRARY})


foreach(LIB IN ITEMS "thermo" "trans")
  add_custom_target(${LIB} ALL
    COMMAND echo ${LIB} | ./${CEA_EXECUTABLE} --legacy-mode > /dev/null
    DEPENDS ${CEA_EXECUTABLE} ${LIB}.inp
  )

  add_custom_command(
    OUTPUT ${LIB}.inp
    COMMAND cp -a ${CMAKE_SOURCE_DIR}/src/${LIB}.inp ${CMAKE_BINARY_DIR}/
  )
endforeach()


install(TARGETS ${CEA_EXECUTABLE} RUNTIME DESTINATION bin)


enable_testing()
add_subdirectory(tests/orig)
add_subdirectory(tests/deton)
add_subdirectory(tests/thermp)
add_subdirectory(tests/shock)
add_subdirectory(tests/rocket)
add_subdirectory(tests/newio)
