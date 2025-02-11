
if(NOT DEFINED CEA_CORE)
  message(FATAL_ERROR "This file must be included from parent CMakeLists.txt.")
endif()


enable_language(CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


set(CEA_OBJECTS cea-objects)

add_library(${CEA_OBJECTS} OBJECT)

target_sources(${CEA_OBJECTS}
  PRIVATE
  cea.cpp
  mod_ffi.f90
)

set_property(TARGET ${CEA_OBJECTS} PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(${CEA_OBJECTS} PRIVATE ${CMAKE_BINARY_DIR})

target_compile_options(${CEA_OBJECTS}
  PRIVATE
  $<$<COMPILE_LANGUAGE:Fortran>:-std=f2023>
)


add_library(${CEA_STATIC_LIBRARY} STATIC $<TARGET_OBJECTS:${CEA_OBJECTS}> $<TARGET_OBJECTS:${CEA_CORE}>)
add_library(${CEA_SHARED_LIBRARY} SHARED $<TARGET_OBJECTS:${CEA_OBJECTS}>)

set_target_properties(${CEA_STATIC_LIBRARY} PROPERTIES OUTPUT_NAME cea)
set_target_properties(${CEA_SHARED_LIBRARY} PROPERTIES OUTPUT_NAME cea)

target_link_libraries(${CEA_STATIC_LIBRARY})
target_link_libraries(${CEA_SHARED_LIBRARY} PRIVATE ${CEA_CORE})

set_target_properties(${CEA_STATIC_LIBRARY} PROPERTIES PUBLIC_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/cea.h)
set_target_properties(${CEA_SHARED_LIBRARY} PROPERTIES PUBLIC_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/cea.h)
