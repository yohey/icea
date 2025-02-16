cmake_minimum_required(VERSION 3.12)

project(CEA-CXX-Examples CXX)

set(EXECUTABLE run_examples)


find_package(CEA REQUIRED)


add_executable(${EXECUTABLE})

target_sources(${EXECUTABLE}
  PRIVATE
  src/main.cpp
  src/example-08.cpp
  src/example-09.cpp
  src/example-10.cpp
  src/example-11.cpp
  src/example-12.cpp
  src/example-13.cpp
)


# target_link_libraries(${EXECUTABLE} PRIVATE cea-static)
target_link_libraries(${EXECUTABLE} PRIVATE cea-shared)
