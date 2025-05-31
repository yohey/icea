cmake_minimum_required(VERSION 3.12)

project(CEA-Fortran-Examples Fortran)

set(EXECUTABLE run_examples)


find_package(CEA REQUIRED)


add_executable(${EXECUTABLE})

target_sources(${EXECUTABLE}
  PRIVATE
  src/main.f90
  src/example-01.f90
  src/example-02.f90
  src/example-03.f90
  src/example-04.f90
  src/example-05.f90
  src/example-06.f90
  src/example-07.f90
  src/example-08.f90
  src/example-09.f90
  src/example-10.f90
  src/example-11.f90
  src/example-12.f90
  src/example-13.f90
)


target_link_libraries(${EXECUTABLE} PRIVATE cea-static)
# target_link_libraries(${EXECUTABLE} PRIVATE cea-shared)
