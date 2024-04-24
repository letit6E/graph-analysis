# Parent project does not export its library target, so this CML implicitly
# depends on being added from it, i.e. the testing is done only from the build
# tree and is not feasible from an install location

project(graph-analysisTests LANGUAGES CXX)

# ---- Dependencies ----

find_package(Catch2 REQUIRED)
include(Catch)

# ---- Tests ----

add_executable(graph-analysis_test source/graph-analysis_test.cpp)
target_link_libraries(
    graph-analysis_test PRIVATE
    graph-analysis_lib
    Catch2::Catch2WithMain
)
target_compile_features(graph-analysis_test PRIVATE cxx_std_20)

catch_discover_tests(graph-analysis_test)

# ---- End-of-file commands ----

add_folders(Test)
