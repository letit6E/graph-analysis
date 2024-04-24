cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_COMPILER "/usr/bin/g++")

include(cmake/prelude.cmake)

project(
    graph-analysis
    VERSION 0.1.0
    DESCRIPTION "Experiments comparing graph analysis algorithms"
    HOMEPAGE_URL "https://github.com/letit6E/graph-analysis"
    LANGUAGES CXX
)

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

# ---- Declare library ----

add_library(
    graph-analysis_lib OBJECT
    source/lib.cpp
)

target_include_directories(
    graph-analysis_lib ${warning_guard}
    PUBLIC
    "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>"
)

target_compile_features(graph-analysis_lib PUBLIC cxx_std_20)

find_package(fmt REQUIRED)
target_link_libraries(graph-analysis_lib PRIVATE fmt::fmt)

# ---- Declare executable ----

add_executable(graph-analysis_exe source/main.cpp)
add_executable(graph-analysis::exe ALIAS graph-analysis_exe)

set_property(TARGET graph-analysis_exe PROPERTY OUTPUT_NAME graph-analysis)

target_compile_features(graph-analysis_exe PRIVATE cxx_std_20)

target_link_libraries(graph-analysis_exe PRIVATE graph-analysis_lib)

# ---- Install rules ----

if(NOT CMAKE_SKIP_INSTALL_RULES)
  include(cmake/install-rules.cmake)
endif()

# ---- Developer mode ----

if(NOT graph-analysis_DEVELOPER_MODE)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  message(
      AUTHOR_WARNING
      "Developer mode is intended for developers of graph-analysis"
  )
endif()

include(cmake/dev-mode.cmake)
