
if(NOT DEFINED CEA_CORE)
  message(FATAL_ERROR "This file must be included from parent CMakeLists.txt.")
endif()


enable_language(CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CEA_MODELICA_LIBRARY cea_modelica)

add_library(${CEA_MODELICA_LIBRARY} STATIC
  $<TARGET_OBJECTS:${CEA_OBJECTS}>
  $<TARGET_OBJECTS:${CEA_CORE}>
  ffi_cea_modelica.cpp
)


if (DEFINED MYCLASSES_MODELICA_INCLUDE)
  if (EXISTS "${MYCLASSES_MODELICA_INCLUDE}/omc/c/ModelicaUtilities.h")
    set(MYCLASSES_MODELICA_C_INCLUDE "${MYCLASSES_MODELICA_INCLUDE}/omc/c")
  else()
    message(FATAL_ERROR "\"${MYCLASSES_MODELICA_INCLUDE}/omc/c/ModelicaUtilities.h\" is not found.")
  endif()
else()
  message(FATAL_ERROR "\"MYCLASSES_MODELICA_INCLUDE\" must be set.")
endif()

target_include_directories(${CEA_MODELICA_LIBRARY}
  PRIVATE
  $<TARGET_PROPERTY:${CEA_STATIC_LIBRARY},SOURCE_DIR>
  ${MYCLASSES_MODELICA_C_INCLUDE}
)


install(TARGETS   ${CEA_MODELICA_LIBRARY}                ARCHIVE DESTINATION lib/omlibrary/CEA/Resources/Library)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/CEA                DESTINATION lib/omlibrary)
install(FILES     ${CMAKE_CURRENT_SOURCE_DIR}/ffi_cea_modelica.h DESTINATION lib/omlibrary/CEA/Resources/Include)
