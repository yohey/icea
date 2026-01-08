
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

target_include_directories(${CEA_MODELICA_LIBRARY} PRIVATE $<TARGET_PROPERTY:${CEA_STATIC_LIBRARY},SOURCE_DIR>)

install(TARGETS   ${CEA_MODELICA_LIBRARY}                ARCHIVE DESTINATION lib/omlibrary/CEA/Resources/Library)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/CEA                DESTINATION lib/omlibrary)
install(FILES     ${CMAKE_CURRENT_SOURCE_DIR}/ffi_cea_modelica.h DESTINATION lib/omlibrary/CEA/Resources/Include)
