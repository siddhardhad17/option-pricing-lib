cmake_minimum_required(VERSION 3.16)

project(optpricing
  VERSION 1.0.0
  DESCRIPTION "C++17 option pricing and SVI volatility surface library"
  LANGUAGES CXX)

# Must come before any target uses CMAKE_INSTALL_INCLUDEDIR in an INSTALL_INTERFACE.
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(OPTPRICING_BUILD_TESTS "Build the unit test suite" ON)
option(OPTPRICING_BUILD_APPS "Build the demo and benchmark applications" ON)
option(OPTPRICING_SANITIZE "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(OPTPRICING_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Threads REQUIRED)

add_library(optpricing
  src/types.cpp
  src/normal.cpp
  src/black_scholes.cpp
  src/barrier.cpp
  src/monte_carlo.cpp
  src/greeks.cpp
  src/nelder_mead.cpp
  src/svi.cpp
  src/no_arbitrage.cpp)
add_library(optpricing::optpricing ALIAS optpricing)

target_include_directories(optpricing
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(optpricing PUBLIC Threads::Threads)
target_compile_features(optpricing PUBLIC cxx_std_17)

# Shared warning and sanitizer policy for every target in this project.
function(optpricing_configure_target target)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wshadow)
    if(OPTPRICING_WARNINGS_AS_ERRORS)
      target_compile_options(${target} PRIVATE -Werror)
    endif()
    if(OPTPRICING_SANITIZE)
      target_compile_options(${target} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
      target_link_options(${target} PRIVATE -fsanitize=address,undefined)
    endif()
  endif()
endfunction()

optpricing_configure_target(optpricing)

if(OPTPRICING_BUILD_TESTS)
  enable_testing()
  add_executable(optpricing_tests
    tests/test_main.cpp
    tests/test_black_scholes.cpp
    tests/test_greeks.cpp
    tests/test_barrier.cpp
    tests/test_monte_carlo.cpp
    tests/test_svi.cpp)
  target_link_libraries(optpricing_tests PRIVATE optpricing)
  optpricing_configure_target(optpricing_tests)
  add_test(NAME optpricing_tests COMMAND optpricing_tests)
endif()

if(OPTPRICING_BUILD_APPS)
  foreach(app price_demo benchmark svi_demo)
    add_executable(${app} apps/${app}.cpp)
    target_link_libraries(${app} PRIVATE optpricing)
    optpricing_configure_target(${app})
  endforeach()
endif()

# ---- Installation and CMake package export -------------------------------------------------
install(TARGETS optpricing
  EXPORT optpricingTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT optpricingTargets
  FILE optpricingTargets.cmake
  NAMESPACE optpricing::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/optpricing)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/optpricingConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/optpricingConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/optpricing)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/optpricingConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/optpricingConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/optpricingConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/optpricing)
