cmake_minimum_required(VERSION 3.18)

project(CEA Fortran)


set(CEA_CORE cea-core)
set(CEA_EXECUTABLE cea)
set(CEA_OBJECTS cea-objects)

set(CEA_Fortran_FLAGS
  -cpp
  -std=f2018
  -DCEA_ROOT='${CMAKE_INSTALL_PREFIX}'
  $<$<CONFIG:Release>:-DNDEBUG>
)

option(WITH_MODELICA "Build with Modelica interfaces" OFF)


add_library(${CEA_CORE} OBJECT)

target_sources(${CEA_CORE}
  PRIVATE
  src/core/mod_cea.f90
  src/core/mod_types.f90
  src/core/mod_constants.f90
  src/core/mod_general.f90
  src/core/mod_functions.f90
  src/core/mod_io.f90
  src/core/mod_legacy_io.f90
  src/interfaces/fortran/mod_ext.f90
  src/interfaces/fortran/cea.f90
)

set_property(TARGET ${CEA_CORE} PROPERTY POSITION_INDEPENDENT_CODE ON)

target_compile_options(${CEA_CORE}
  PRIVATE
  ${CEA_Fortran_FLAGS}
)


add_executable(${CEA_EXECUTABLE})

target_sources(${CEA_EXECUTABLE}
  PRIVATE
  src/main.f90
)

target_compile_options(${CEA_EXECUTABLE}
  PRIVATE
  ${CEA_Fortran_FLAGS}
)

target_link_libraries(${CEA_EXECUTABLE} PRIVATE ${CEA_CORE})


if(WIN32)
  set(PATH_DELIMITER "\\")
  set(NULL_DEVICE "NUL")
else()
  set(PATH_DELIMITER "/")
  set(NULL_DEVICE "/dev/null")
endif()

foreach(LIB IN ITEMS "thermo" "trans")
  add_custom_target(${LIB}.lib ALL
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/src/core/${LIB}.inp ${CMAKE_BINARY_DIR}/
    COMMAND ${CMAKE_COMMAND} -E remove -f ${LIB}.lib
    COMMAND echo ${LIB} | .${PATH_DELIMITER}${CEA_EXECUTABLE} --legacy-mode > ${NULL_DEVICE}
    DEPENDS ${CEA_EXECUTABLE}
  )
endforeach()


set(CEA_STATIC_LIBRARY cea-static)
set(CEA_SHARED_LIBRARY cea-shared)

add_subdirectory(src/interfaces/c++)

add_subdirectory(src/interfaces/python)

if(WITH_MODELICA)
  add_subdirectory(src/interfaces/modelica)
endif()


install(TARGETS ${CEA_EXECUTABLE})

install(TARGETS ${CEA_STATIC_LIBRARY} ${CEA_SHARED_LIBRARY}
  EXPORT                    cea-config
  LIBRARY DESTINATION       lib
  INCLUDES DESTINATION      include include/cea
  PUBLIC_HEADER DESTINATION include/cea
)

install(FILES ${CMAKE_BINARY_DIR}/thermo.lib DESTINATION lib/cea/)
install(FILES ${CMAKE_BINARY_DIR}/trans.lib  DESTINATION lib/cea/)
install(FILES ${CMAKE_BINARY_DIR}/cea.mod    DESTINATION include/cea/)

install(FILES tests/orig/cea2.inp  DESTINATION share/doc/cea/examples/orig/)
install(FILES src/core/thermo.inp  DESTINATION share/doc/cea/examples/orig/)
install(FILES src/core/trans.inp   DESTINATION share/doc/cea/examples/orig/)


file(GLOB EXAMPLES_FORTRAN examples/fortran/src/*.f90)
file(GLOB EXAMPLES_CXX     examples/c++/src/*.cpp)
file(GLOB EXAMPLES_PYTHON  examples/python/*.py)

install(FILES     examples/fortran/CMakeLists.txt DESTINATION share/doc/cea/examples/fortran/)
install(FILES     examples/c++/CMakeLists.txt DESTINATION share/doc/cea/examples/c++/)
install(FILES     ${EXAMPLES_FORTRAN} DESTINATION share/doc/cea/examples/fortran/src/)
install(FILES     ${EXAMPLES_CXX}     DESTINATION share/doc/cea/examples/c++/src/)
install(FILES     ${EXAMPLES_PYTHON}  DESTINATION share/doc/cea/examples/python/)

install(EXPORT cea-config DESTINATION share/cmake/cea/ EXPORT_LINK_INTERFACE_LIBRARIES)

add_subdirectory(src/docs)


enable_testing()
add_subdirectory(tests/orig)
add_subdirectory(tests/deton)
add_subdirectory(tests/thermp)
add_subdirectory(tests/shock)
add_subdirectory(tests/rocket)
add_subdirectory(tests/interfaces/fortran)
add_subdirectory(tests/interfaces/c++)
add_subdirectory(tests/interfaces/python)
