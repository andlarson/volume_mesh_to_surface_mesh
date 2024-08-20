cmake_minimum_required(VERSION 3.28)
project(volume_mesh_to_surface_mesh LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# ------------------------------------------------------------------------------
#                                   Executable 
# ------------------------------------------------------------------------------

set(CMAKE_MODULE_PATH "/Users/andrewlarson/Desktop/misc/programming/cmake_utils/find_modules")

find_package(CLI11 REQUIRED)
find_package(Gmsh MODULE QUIET REQUIRED)

add_executable(volume_mesh_to_surface_mesh "${CMAKE_SOURCE_DIR}/src/main.cpp")

target_link_libraries(volume_mesh_to_surface_mesh
                      PUBLIC
                      Gmsh::Gmsh
                      CLI11::CLI11)
