cmake_minimum_required(VERSION 3.12)

project(CEA Fortran)


set(CEA_CORE cea-core)
set(CEA_EXECUTABLE cea)


add_library(${CEA_CORE} OBJECT)

target_sources(${CEA_CORE}
  PRIVATE
  src/core/mod_cea_core.f90
  src/core/mod_cea_types.f90
  src/core/mod_constants.f90
  src/core/mod_general.f90
  src/core/mod_functions.f90
  src/core/mod_io.f90
  src/core/mod_legacy_io.f90
  src/interfaces/fortran/cea.f90
)

set_property(TARGET ${CEA_CORE} PROPERTY POSITION_INDEPENDENT_CODE ON)

target_compile_options(${CEA_CORE}
  PRIVATE
  "-cpp"
  "-DCEA_ROOT='${CMAKE_INSTALL_PREFIX}'"
  "-std=f2018"
)


add_executable(${CEA_EXECUTABLE})

target_sources(${CEA_EXECUTABLE}
  PRIVATE
  src/main.f90
)

target_compile_options(${CEA_EXECUTABLE}
  PRIVATE
  "-cpp"
  "-std=f2018"
)

target_link_libraries(${CEA_EXECUTABLE} PRIVATE ${CEA_CORE})


foreach(LIB IN ITEMS "thermo" "trans")
  add_custom_target(${LIB}.lib ALL
    COMMAND cp -a ${CMAKE_SOURCE_DIR}/src/core/${LIB}.inp ${CMAKE_BINARY_DIR}/
    COMMAND rm -f ${LIB}.lib
    COMMAND echo ${LIB} | ./${CEA_EXECUTABLE} --legacy-mode > /dev/null
    DEPENDS ${CEA_EXECUTABLE}
  )
endforeach()


set(CEA_STATIC_LIBRARY cea-static)
set(CEA_SHARED_LIBRARY cea-shared)

add_subdirectory(src/interfaces/c++)


install(TARGETS ${CEA_EXECUTABLE})

install(TARGETS ${CEA_STATIC_LIBRARY}
  EXPORT               cea-config
  LIBRARY DESTINATION  lib
  INCLUDES DESTINATION include
)

install(TARGETS ${CEA_SHARED_LIBRARY}
  EXPORT               cea-config
  LIBRARY DESTINATION  lib
  INCLUDES DESTINATION include
)

install(FILES ${CMAKE_BINARY_DIR}/thermo.lib DESTINATION lib/)
install(FILES ${CMAKE_BINARY_DIR}/trans.lib  DESTINATION lib/)
install(FILES tests/orig/cea2.inp DESTINATION share/examples/cea/orig/)
install(FILES src/core/thermo.inp DESTINATION share/examples/cea/orig/)
install(FILES src/core/trans.inp  DESTINATION share/examples/cea/orig/)

install(FILES ${CMAKE_BINARY_DIR}/cea.mod DESTINATION include/)

install(EXPORT cea-config DESTINATION share/cmake/cea/ EXPORT_LINK_INTERFACE_LIBRARIES)


enable_testing()
add_subdirectory(tests/orig)
add_subdirectory(tests/deton)
add_subdirectory(tests/thermp)
add_subdirectory(tests/shock)
add_subdirectory(tests/rocket)
add_subdirectory(tests/interfaces/fortran)
# add_subdirectory(tests/c++)
