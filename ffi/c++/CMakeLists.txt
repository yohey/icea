
enable_language(CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


add_library(${CEA_CXX_LIBRARY})

set_target_properties(${CEA_CXX_LIBRARY}
  PROPERTIES
  OUTPUT_NAME cea
)

target_sources(${CEA_CXX_LIBRARY}
  PRIVATE
  src/cea.cpp
  src/mod_ffi.f90
)

target_include_directories(${CEA_CXX_LIBRARY} PRIVATE ${CMAKE_BINARY_DIR})

target_link_libraries(${CEA_CXX_LIBRARY} PRIVATE ${CEA_CORE})

target_compile_options(${CEA_CXX_LIBRARY}
  PRIVATE
  $<$<COMPILE_LANGUAGE:Fortran>:-std=f2023>
)

