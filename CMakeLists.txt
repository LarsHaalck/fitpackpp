cmake_minimum_required(VERSION 3.5)
project(fitpackpp VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)
option(BUILD_EXAMPLES "Build examples" OFF)

add_subdirectory(fitpack)
add_subdirectory(fitpackpp)

if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

# install targets
include(GNUInstallDirs)
install(
  EXPORT fitpackpp-targets
  FILE fitpackppTargets.cmake
  NAMESPACE fitpackpp::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fitpackpp
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/fitpackppConfig.cmake"
  INSTALL_DESTINATION lib/cmake/fitpackpp
)

install(
  FILES
  "${CMAKE_CURRENT_BINARY_DIR}/fitpackppConfig.cmake"
  DESTINATION lib/cmake/fitpackpp
)
