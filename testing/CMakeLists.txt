# Cmake file for testing suite using googletest
# Setup from https://google.github.io/googletest/quickstart-cmake.html

# Required by googletest
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Pull googletest dependencies, uses GIT_* tags instead of URL, as url didn't pull properly
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest
  GIT_TAG b514bdc898e2951020cbdca1304b75f5950d1f59 # refs/tags/v1.15.2
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

# Testing executable, update for more tests and dependencies 
add_executable(
  ${PROJECT_NAME}_testing

  # Test files
  test_settings.cpp
  test_discretization.cpp
  test_solver.cpp
  test_storage.cpp

  # Dependecies
  ../src/settings/settings.cpp
)

# Include source directory to be visible from everytest
target_include_directories(${PROJECT_NAME}_testing PUBLIC ${PROJECT_SOURCE_DIR}/src)

# Link external library googletest
target_link_libraries(
  ${PROJECT_NAME}_testing
  GTest::gtest_main
)

# Discover google tests
include(GoogleTest)
gtest_discover_tests(${PROJECT_NAME}_testing)