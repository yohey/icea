cmake_minimum_required(VERSION 3.12)

project(CEA-CXX-Examples CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


set(EXECUTABLE run_examples)


find_package(CEA REQUIRED)


add_executable(${EXECUTABLE})

target_sources(${EXECUTABLE}
  PRIVATE
  src/main.cpp
  src/example-01.cpp
  src/example-02.cpp
  src/example-03.cpp
  src/example-04.cpp
  src/example-05.cpp
  src/example-06.cpp
  src/example-07.cpp
  src/example-08.cpp
  src/example-09.cpp
  src/example-10.cpp
  src/example-11.cpp
  src/example-12.cpp
  src/example-13.cpp
  src/example-14.cpp
)


target_link_libraries(${EXECUTABLE} PRIVATE cea-shared)
