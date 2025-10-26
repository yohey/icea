
if(NOT DEFINED CEA_SHARED_LIBRARY)
  message(FATAL_ERROR "This file must be included from parent CMakeLists.txt.")
endif()

set(SPHINX_BUILD_EXECUTABLE  "sphinx-build"  CACHE STRING "Path to sphinx-build executable")

set(SPHINX_ORIG_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_COPY_SOURCE_DIR ${CMAKE_BINARY_DIR}/docs/source)
set(SPHINX_HTML_OUTPUT_DIR ${CMAKE_BINARY_DIR}/docs/html)
set(PYTHON_PACKAGE_DIR ${CMAKE_SOURCE_DIR}/src/interfaces/python/icea)

add_custom_command(
  OUTPUT ${SPHINX_COPY_SOURCE_DIR}/index.rst
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${SPHINX_ORIG_SOURCE_DIR} ${SPHINX_COPY_SOURCE_DIR}
  COMMENT "Copying Sphinx documentation source to build directory"
  VERBATIM
)

add_custom_target(copy_sphinx_sources
  DEPENDS ${SPHINX_COPY_SOURCE_DIR}/index.rst
)

set(SPHINX_ENV_DEFINITION "PYTHONPATH=${CMAKE_SOURCE_DIR}/src/interfaces/python")

if(APPLE)
  set(SPHINX_ENV_DEFINITION "${SPHINX_ENV_DEFINITION};DYLD_LIBRARY_PATH=${CMAKE_BINARY_DIR}")
elseif(UNIX)
  set(SPHINX_ENV_DEFINITION "${SPHINX_ENV_DEFINITION};LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}")
elseif(WIN32)
  set(SPHINX_ENV_DEFINITION "${SPHINX_ENV_DEFINITION};PATH=${CMAKE_BINARY_DIR}\;$ENV{PATH}")
endif()

add_custom_target(doc
  COMMAND ${CMAKE_COMMAND} -E env
          ${SPHINX_ENV_DEFINITION}
          ${SPHINX_BUILD_EXECUTABLE} -q
            -b html
            ${SPHINX_COPY_SOURCE_DIR}
            ${SPHINX_HTML_OUTPUT_DIR}
  DEPENDS copy_sphinx_sources ${CEA_SHARED_LIBRARY}
  COMMENT "Building HTML documentation with sphinx-build"
  VERBATIM
)

set_directory_properties(PROPERTIES
    ADDITIONAL_MAKE_CLEAN_FILES "${SPHINX_HTML_OUTPUT_DIR}"
)
