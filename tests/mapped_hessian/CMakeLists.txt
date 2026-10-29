cmake_minimum_required(VERSION 3.15)
project(fusionio_mapped_hessian_regression LANGUAGES C CXX)
enable_testing()

# Standalone regression: no changes to the main build, installation, or defaults.
get_filename_component(_default_source "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
set(FUSIONIO_SOURCE_DIR "${_default_source}" CACHE PATH "FusionIO source tree under test")
find_package(HDF5 REQUIRED COMPONENTS C)
add_executable(mapped_hessian_regression
  mapped_hessian_regression.cpp
  "${FUSIONIO_SOURCE_DIR}/m3dc1_lib/m3dc1_mesh.cpp"
  "${FUSIONIO_SOURCE_DIR}/m3dc1_lib/m3dc1_field.cpp"
  "${FUSIONIO_SOURCE_DIR}/m3dc1_lib/m3dc1_coord_map.cpp"
  "${FUSIONIO_SOURCE_DIR}/m3dc1_lib/m3dc1_stell_field.cpp")
target_compile_features(mapped_hessian_regression PRIVATE cxx_std_17)
target_include_directories(mapped_hessian_regression PRIVATE
  "${FUSIONIO_SOURCE_DIR}/m3dc1_lib"
  "${FUSIONIO_SOURCE_DIR}/output_stream"
  ${HDF5_INCLUDE_DIRS})
target_link_libraries(mapped_hessian_regression PRIVATE ${HDF5_LIBRARIES})
add_test(NAME mapped_hessian_regression COMMAND mapped_hessian_regression)

